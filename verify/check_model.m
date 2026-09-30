function ok = check_model(root)
%CHECK_MODEL  Check the signals and the recursion against the paper and exact results.
%   ok = check_model(root) prints one PASS/FAIL line per check. None of these checks
%   uses the Python code: the references are the paper and closed-form results.

fprintf('\n== Model checks (against the paper and exact results) ==\n');
rng(7);
ok = true;

% --- h: the impulse response of Fig. 2 ---------------------------------------
h = load(fullfile(root, 'data', 'h.csv'), '-ascii');
[pk, ipk] = max(h);
ok = check('h: M = 128 taps, Eq. (88)', numel(h), numel(h) == 128) && ok;
ok = check('h: peak value, Fig. 2 reads 0.050', pk, abs(pk - 0.050) < 0.001) && ok;
ok = check('h: peak at tap 33 counting from 0, Fig. 2', ipk - 1, ipk - 1 == 33) && ok;

% --- AR(1) input, Sec. 4 -----------------------------------------------------
a = 0.9;
v_u = 1;
x = af.ar1_input(200000, 1, a, v_u);
var_x = v_u / (1 - a ^ 2);
ok = check('input: var(x) / (v_u / (1 - a^2)) = 1', var(x) / var_x, ...
           abs(var(x) / var_x - 1) < 0.03) && ok;
r1 = (x(1:end-1)' * x(2:end)) / (x' * x);
ok = check('input: lag-1 correlation = -a = -0.9', r1, abs(r1 + a) < 0.01) && ok;

% --- output power for the SNR, Eq. (89) --------------------------------------
xh = filter(h, 1, x);                           % x_t' h for every t
p_emp = mean(xh(numel(h):end) .^ 2);
p_exact = af.noise_var_for_snr(h, 0, a, v_u);   % SNR = 0 dB gives h' R h
ok = check('SNR: sampled E[(x''h)^2] / (h'' R h) = 1', p_emp / p_exact, ...
           abs(p_emp / p_exact - 1) < 0.03) && ok;

% --- generalized Gaussian noise, Eqs. (59)-(62) ------------------------------
% Kolmogorov-Smirnov distance between the sampled and the exact CDF of |eta|:
% P(|eta| <= z) = P(G <= (z/c)^beta) with G ~ Gamma(1/beta, 1), c = kappa(beta)
% for unit variance. The threshold is the 1% critical value.
n = 100000;
for beta_star = [0.2 1 2]
  e = af.gg_rnd([n 1], 1, beta_star);
  c = af.gg_kappa(beta_star);
  z = sort(abs(e));
  F = gammainc((z / c) .^ beta_star, 1 / beta_star);
  D = max(max((1:n)' / n - F), max(F - (0:n-1)' / n));
  ok = check(sprintf('noise: KS distance to exact CDF, beta* = %.1f', beta_star), D, ...
             D < 1.63 / sqrt(n)) && ok;
  if beta_star >= 1   % at beta* = 0.2 the kurtosis is ~770 and var() is too noisy
    ok = check(sprintf('noise: sample variance = 1, beta* = %.1f', beta_star), var(e), ...
               abs(var(e) - 1) < 0.03) && ok;
  end
end

% --- tau, Eq. (69) -----------------------------------------------------------
sigma = 0.3;
ok = check('tau(beta = 2) / v_eta = 1, Sec. 2.5', af.gg_tau(sigma, 2) / sigma ^ 2, ...
           abs(af.gg_tau(sigma, 2) / sigma ^ 2 - 1) < 1e-12) && ok;
ok = check('tau(beta = 1) / (sigma / sqrt(2)) = 1, not Eq. (80)', ...
           af.gg_tau(sigma, 1) / (sigma / sqrt(2)), ...
           abs(af.gg_tau(sigma, 1) / (sigma / sqrt(2)) - 1) < 1e-12) && ok;

% --- KF, beta = 2, eps = 0 is the exact Bayesian posterior, Sec. 2.5 --------
% With Gaussian noise, a constant state and V_0 = v0 I, the recursion (26)-(33)
% must return the posterior mean (X'X + v_eta/v0 I)^-1 X'y, whatever the data.
T = 400;
M = numel(h);
v_eta = af.noise_var_for_snr(h, 5, a, v_u);
x = af.ar1_input(T, 1, a, v_u);
eta = sqrt(v_eta) * randn(T, 1);
opt = struct('eps', 0, 'sigma_eta', sqrt(v_eta));
out = af.run_filter(h, x, eta, 'KF', 2, opt);
Xmat = toeplitz(x, [x(1) zeros(1, M - 1)]);     % row t is x_t', zeros before t = 1
y = Xmat * h + eta;
v0 = (h' * h) / M;
w_exact = (Xmat' * Xmat + v_eta / v0 * eye(M)) \ (Xmat' * y);
err = max(abs(out.w - w_exact)) / max(abs(w_exact));
ok = check('KF, beta = 2, eps = 0 vs closed-form posterior', err, err < 1e-8) && ok;

if ok
  fprintf('  all model checks passed\n');
else
  fprintf('  SOME MODEL CHECKS FAILED\n');
end
end

function pass = check(label, value, pass)
if pass
  s = 'PASS';
else
  s = 'FAIL';
end
fprintf('  %-56s %12.5g   %s\n', label, value, s);
end
