function k = paper_params(beta, target_db, n_iter)
%PAPER_PARAMS  Knob values published in Tables 1 and 2 of the paper.
%   k = af.paper_params(beta, target_db, n_iter) returns a struct with
%     k.SG.mu, k.fKF.reg (= tau/v), k.sKF.eps, k.vKF.eps, k.KF.eps
%
%   Table 1: beta = 2 and beta = 1, target -20 dB, I = 0 (Fig. 3).
%   Table 2: beta = 1, targets -15, -20, -25 dB, I = 0 and 1 (Fig. 4).
%   vKF is not in the tables (the paper leaves it out of its figures because it
%   matched sKF), so it gets the sKF value. SG has no iterations, so its entry is
%   the same for I = 0 and I = 1.

if beta == 2
  if target_db ~= -20 || n_iter ~= 0
    error('af:paper_params', 'Table 1 gives beta = 2 only for -20 dB and I = 0');
  end
  knobs = [1.1e-4; 8.2e3; 3.2e-10; 3.6e-11];
elseif beta == 1
  % Table 2. Its -20 dB, I = 0 column is also the beta = 1 column of Table 1.
  %        -15 dB           -20 dB           -25 dB
  %        I = 0   I = 1    I = 0   I = 1    I = 0   I = 1
  tab = [5.4e-5  5.4e-5   2.7e-5  2.7e-5   1.4e-5  1.4e-5;    % SG   mu
         3.5e3   5.0e3    1.1e4   1.6e4    3.4e4   4.3e4;     % fKF  tau/v
         1.0e-7  7.1e-8   2.7e-8  2.2e-8   7.7e-9  6.6e-9;    % sKF  eps
         7.3e-8  6.0e-8   2.2e-8  2.2e-8   6.0e-9  6.0e-9];   % KF   eps
  col = find([-15 -20 -25] == target_db);
  if isempty(col) || ~any(n_iter == [0 1])
    error('af:paper_params', 'Table 2 has targets -15, -20, -25 dB and I = 0, 1');
  end
  knobs = tab(:, 2 * (col - 1) + n_iter + 1);
else
  error('af:paper_params', 'the tables only cover beta = 2 and beta = 1');
end

k.SG = struct('mu', knobs(1));
k.fKF = struct('reg', knobs(2));
k.sKF = struct('eps', knobs(3));
k.vKF = struct('eps', knobs(3));
k.KF = struct('eps', knobs(4));
end
