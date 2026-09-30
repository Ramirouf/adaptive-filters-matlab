% RUN_ALL  Checks, cross-validation, Figs. 3-4 of Szczecinski, Benesty & Kuhn (2026), and a
%          comparison with the published Fig. 3 curves.
%
%   octave --no-gui run_all.m        full run, N = 100 as in the paper
%   >> QUICK = true; run_all          N = 10, fits MATLAB Online's 15-minute compute limit
%
%   Writes figures/*.png, results/fig*.csv and results/summary.csv.

if ~exist('QUICK', 'var')
  QUICK = false;
end
if QUICK
  N = 10;
else
  N = 100;
end
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'sim'), fullfile(root, 'verify'));
t_start = tic;

% ---- checks first: seconds, and they catch most mistakes ---------------------
ok_model = check_model(root);
ok_python = compare_python(root);

% ---- the experiment of Sec. 4 -----------------------------------------------
setup.h = load(fullfile(root, 'data', 'h.csv'), '-ascii');
setup.a = 0.9;             % AR(1) coefficient
setup.v_u = 1;             % variance of u_t; not stated in the paper
setup.snr_db = 5;          % Eq. (89)
setup.beta_star = 0.2;     % shape of the actual noise
setup.fs = 8000;           % Hz
setup.seed = 1;

fig3 = {'SG', 0; 'fKF', 0; 'sKF', 0; 'KF', 0};
fig4 = {'SG', 0; 'fKF', 0; 'fKF', 1; 'sKF', 0; 'sKF', 1; 'KF', 0; 'KF', 1};

fprintf('\n== Figs. 3 and 4, N = %d ==\n', N);
stats = [make_figure('fig3a', 2, -20, fig3, 12.0, N, setup, root), ...
         make_figure('fig3b', 1, -20, fig3, 1.2, N, setup, root), ...
         make_figure('fig4a', 1, -15, fig4, 1.0, N, setup, root), ...
         make_figure('fig4b', 1, -20, fig4, 1.2, N, setup, root), ...
         make_figure('fig4c', 1, -25, fig4, 2.0, N, setup, root)];

fid = fopen(fullfile(root, 'results', 'summary.csv'), 'w');
fprintf(fid, '%% N = %d; floor = mean of the last 10%% (dB); t_conv = time to within 1 dB of the floor (s)\n', N);
fprintf(fid, '%% panel, variant, beta, I, target_db, floor_db, t_conv_s, run_s\n');
for i = 1:numel(stats)
  s = stats(i);
  fprintf(fid, '%s, %s, %d, %d, %d, %.2f, %.4f, %.1f\n', s.panel, s.variant, s.beta, ...
          s.n_iter, s.target_db, s.floor_db, s.t_conv, s.run_s);
end
fclose(fid);

% ---- the paper's headline: robust filters converge ~10x faster (Fig. 3) -----
fprintf('\n== Speed-up of beta = 1 over beta = 2 at the -20 dB target (Fig. 3) ==\n');
fprintf('  %-5s %14s %14s %14s %14s %9s\n', 'filter', 'floor b=2 [dB]', 'floor b=1 [dB]', ...
        't_conv b=2 [s]', 't_conv b=1 [s]', 'speed-up');
variants = {'SG', 'fKF', 'sKF', 'KF'};
for i = 1:numel(variants)
  a2 = stats(strcmp({stats.panel}, 'fig3a') & strcmp({stats.variant}, variants{i}));
  a1 = stats(strcmp({stats.panel}, 'fig3b') & strcmp({stats.variant}, variants{i}));
  fprintf('  %-5s %14.2f %14.2f %14.3f %14.3f %8.1fx\n', variants{i}, a2.floor_db, ...
          a1.floor_db, a2.t_conv, a1.t_conv, a2.t_conv / a1.t_conv);
end

% ---- against the published curves: markers read off Figs. 3a and 3b ---------
ok_paper = compare_paper(root);

fprintf('\nmodel checks: %s   cross-validation: %s   paper markers: %s   total time %.0f s\n', ...
        mat2str(ok_model), mat2str(ok_python), mat2str(ok_paper), toc(t_start));
