function out = run_filter(h, x, eta, variant, beta, opt)
%RUN_FILTER  The unified adaptive filter of Secs. 2.2-2.4 and 3.1, over N realizations.
%   out = af.run_filter(h, x, eta, variant, beta, opt)
%
%   h        M x 1   impulse response to identify, Eq. (88)
%   x        T x N   input signal, one realization per column
%   eta      T x N   measurement noise, same layout; y_t = x_t' h + eta_t, Eq. (86)
%   variant  'SG', 'fKF', 'sKF', 'vKF' or 'KF'
%   beta     noise shape the filter ASSUMES: 2 = conventional, 1 = robust
%   opt      struct with fields
%     mu            SG step size, Eqs. (54), (65)
%     reg           fKF regularization tau/v, Eqs. (56), (68)
%     eps           sKF/vKF/KF state-noise variance, Eqs. (26), (44)
%     sigma_eta     noise std. dev. the filter assumes (sKF/vKF/KF only)
%     v0            prior variance, V_0 = v0 I. Default ||h||^2 / M. The paper
%                   (p. 3) leaves it as a free parameter and gives no value.
%     n_iter        I, extra iterations of the gain, Eqs. (37)-(40). Default 0.
%     record_every  store the misalignment every k samples. Default 1.
%     w_every       also store the weights every k samples. Default 0 (off).
%
%   out.mis     misalignment ||w_t - h||^2 / ||h||^2, Eq. (90), averaged over the
%               N realizations, Eq. (91). Linear scale.
%   out.t       sample index of each stored value, counting from 0
%   out.w       M x N final estimates
%   out.w_traj  M x N x K weights at samples 0, w_every, 2 w_every, ... ([] if off)
%
%   All five variants are one recursion. With the prediction error
%       e_t     = y_t - x_t' w_{t-1}                                (29)
%   the gain multiplier and the update are
%       alpha_t = 1 / (tau |e_t|^(2-beta) + s_t)                    (31) with (64), = (75)
%       w_t     = w_{t-1} + kappa_t alpha_t e_t                     (32)
%   and the variants differ only in the covariance model (41), which sets kappa_t
%   and s_t:
%       KF    V_t full       kappa = V x      s = x' kappa          (27), (28)
%       vKF   diag(v_t)      kappa = v .* x   s = x' kappa          (45), (46)
%       sKF   v_t I          kappa = v x      s = v ||x||^2         (45), (46)
%       fKF   v I, frozen    kappa = x        s = ||x||^2           tau -> tau/v, (68)
%       SG    s_t -> 0       w_t = w_{t-1} + mu x |e_t|^(beta-1) sign(e_t)   (65)
%   For fKF, dividing kappa, s and tau by the constant v leaves kappa*alpha
%   unchanged, which is why its only knob is the ratio tau/v.

[T, N] = size(x);
if ~isequal(size(eta), [T N])
  error('af:run_filter', 'x and eta must have the same size');
end
h = h(:);
M = numel(h);
hnorm2 = h' * h;

n_iter = get_opt(opt, 'n_iter', 0);
record_every = get_opt(opt, 'record_every', 1);
w_every = get_opt(opt, 'w_every', 0);

is_sg = strcmp(variant, 'SG');
is_f = strcmp(variant, 'fKF');
is_s = strcmp(variant, 'sKF');
is_v = strcmp(variant, 'vKF');
is_k = strcmp(variant, 'KF');
if is_sg
  mu = opt.mu;
elseif is_f
  d = opt.reg;                                  % tau / v, Eq. (68)
elseif is_s || is_v || is_k
  epsv = opt.eps;
  d = af.gg_tau(opt.sigma_eta, beta);           % Eq. (69)
  v0 = get_opt(opt, 'v0', hnorm2 / M);
else
  error('af:run_filter', 'unknown variant ''%s''', variant);
end

n_rec = ceil(T / record_every);
mis_sum = zeros(n_rec, 1);
w = zeros(M, N);
if w_every > 0
  w_traj = zeros(M, N, ceil(T / w_every));
else
  w_traj = [];
end

% The full-covariance KF runs one realization at a time: in Octave that is faster
% than keeping N covariance matrices in a 3-D array. The others run all N at once.
if is_k
  blocks = num2cell(1:N);
else
  blocks = {1:N};
end

for b = 1:numel(blocks)
  cols = blocks{b};
  B = numel(cols);
  xb = x(:, cols);
  eb = eta(:, cols);
  W = zeros(M, B);                              % w_0 = 0
  X = zeros(M, B);                              % x_t = [x_t, ..., x_{t-M+1}]', Eq. (87)
  if is_s
    v = v0 * ones(1, B);                        % one variance per realization
  elseif is_v
    v = v0 * ones(M, B);                        % one variance per tap
  elseif is_k
    V = v0 * eye(M);                            % V_0 = v0 I (here B = 1)
  end
  rec = 0;
  rec_w = 0;

  for t = 1:T
    X = [xb(t, :); X(1:M-1, :)];
    y = h' * X + eb(t, :);                      % Eq. (86)
    e0 = y - sum(X .* W, 1);                    % Eq. (29)

    if is_sg
      W = W + X .* (mu * abs(e0) .^ (beta - 1) .* sign(e0));   % Eq. (65)
    else
      % prediction and preliminary Kalman gain
      if is_s
        v = v + epsv;                           % Eq. (44)
        kap = X .* v;                           % Eq. (45)
        s = v .* sum(X .^ 2, 1);                % Eq. (46)
      elseif is_v
        v = v + epsv;                           % Eq. (44)
        kap = v .* X;                           % Eq. (45)
        s = sum(X .* kap, 1);                   % Eq. (46)
      elseif is_k
        V(1:M+1:end) = V(1:M+1:end) + epsv;     % Eq. (26)
        kap = V * X;                            % Eq. (27)
        s = X' * kap;                           % Eq. (28)
      else
        kap = X;                                % fKF, Eq. (68)
        s = sum(X .^ 2, 1);
      end

      % gain multiplier and update, refined I times, Eqs. (37)-(40)
      e_i = e0;
      for i = 0:n_iter
        alpha = 1 ./ (d * abs(e_i) .^ (2 - beta) + s);   % Eqs. (31), (75)
        W_new = W + kap .* (alpha .* e0);                 % Eq. (38): e_{t,0} in the numerator
        e_i = y - sum(X .* W_new, 1);
      end
      W = W_new;                                % Eq. (40)

      % posterior variance
      if is_s
        v = v .* (1 - s .* alpha / M);          % Eq. (47)
      elseif is_v
        v = v .* (1 - kap .* X .* alpha);       % Eq. (47)
      elseif is_k
        V = V - (alpha * kap) * kap';           % Eq. (33), using V x = kappa
      end
    end

    if mod(t - 1, record_every) == 0
      rec = rec + 1;
      dw = W - h;
      mis_sum(rec) = mis_sum(rec) + sum(dw(:) .^ 2);
    end
    if w_every > 0 && mod(t - 1, w_every) == 0
      rec_w = rec_w + 1;
      w_traj(:, cols, rec_w) = W;
    end
  end
  w(:, cols) = W;
end

out.mis = mis_sum / (N * hnorm2);
out.t = (0:n_rec-1)' * record_every;
out.w = w;
out.w_traj = w_traj;
end

function val = get_opt(opt, name, default)
if isfield(opt, name) && ~isempty(opt.(name))
  val = opt.(name);
else
  val = default;
end
end
