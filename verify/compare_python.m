function ok = compare_python(root)
%COMPARE_PYTHON  Cross-validate af.run_filter against the Python reproduction.
%   ok = compare_python(root) loads verify/reference/ (written by export_reference.py),
%   runs the same 14 configurations on the same input and noise, and prints a
%   PASS/FAIL table.
%
%   err_w  largest relative gap between the weight vectors, over all stored snapshots:
%          max_k ||w_k(MATLAB) - w_k(Python)||_inf / ||w_k(Python)||_inf
%   err_m  largest relative gap between the misalignment curves, over every sample
%
%   Tolerance 1e-9: rounding differences from a different order of operations stay
%   far below it, while any real discrepancy lands far above it.

tol = 1e-9;
ref = fullfile(root, 'verify', 'reference');
fprintf('\n== Cross-validation against the Python reproduction ==\n');
if ~exist(fullfile(ref, 'x.csv'), 'file')
  fprintf('  no reference in verify/reference: run verify/export_reference.py first\n');
  ok = false;
  return;
end

s = load(fullfile(ref, 'setup.csv'), '-ascii');   % T N w_every a v_u snr_db beta_star v_eta
w_every = s(3);
h = load(fullfile(ref, 'h.csv'), '-ascii');
x = load(fullfile(ref, 'x.csv'), '-ascii');
eta = load(fullfile(ref, 'eta.csv'), '-ascii');
fprintf('  %d samples x %d realizations, weights compared every %d samples\n', ...
        size(x, 1), size(x, 2), w_every);

v_eta = af.noise_var_for_snr(h, s(6), s(4), s(5));
err_v = abs(v_eta - s(8)) / s(8);
ok = err_v < tol;
fprintf('  noise variance for SNR = %g dB, Eq. (89): relative gap %.1e   %s\n', ...
        s(6), err_v, pass_str(err_v < tol));

configs = {'SG', 2, 0; 'fKF', 2, 0; 'sKF', 2, 0; 'vKF', 2, 0; 'KF', 2, 0; ...
           'SG', 1, 0; 'fKF', 1, 0; 'sKF', 1, 0; 'vKF', 1, 0; 'KF', 1, 0; ...
           'fKF', 1, 1; 'sKF', 1, 1; 'vKF', 1, 1; 'KF', 1, 1};
fprintf('\n  %-8s %4s %3s %12s %12s   %s\n', 'variant', 'beta', 'I', 'err_w', 'err_m', 'result');
worst = 0;
for c = 1:size(configs, 1)
  [variant, beta, n_iter] = configs{c, :};
  name = sprintf('%s_beta%d_I%d', variant, beta, n_iter);
  k = af.paper_params(beta, -20, n_iter);
  opt = k.(variant);
  opt.sigma_eta = sqrt(v_eta);
  opt.n_iter = n_iter;
  opt.w_every = w_every;
  out = af.run_filter(h, x, eta, variant, beta, opt);

  mis_py = load(fullfile(ref, [name '_mis.csv']), '-ascii');
  w_py = load(fullfile(ref, [name '_w.csv']), '-ascii');       % K x (N M)
  w_ml = reshape(out.w_traj, [], size(out.w_traj, 3))';        % same layout
  err_w = max(max(abs(w_ml - w_py), [], 2) ./ max(abs(w_py), [], 2));
  err_m = max(abs(out.mis - mis_py) ./ mis_py);

  pass = err_w < tol && err_m < tol;
  ok = ok && pass;
  worst = max([worst err_w err_m]);
  fprintf('  %-8s %4d %3d %12.2e %12.2e   %s\n', variant, beta, n_iter, err_w, err_m, pass_str(pass));
end

fprintf('\n  largest gap over all configurations: %.2e (tolerance %.0e)\n', worst, tol);
if ok
  fprintf('  cross-validation PASSED\n');
else
  fprintf('  cross-validation FAILED\n');
end
end

function s = pass_str(pass)
if pass
  s = 'PASS';
else
  s = 'FAIL';
end
end
