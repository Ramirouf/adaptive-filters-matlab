function stats = make_figure(name, beta, target_db, curves, dur_s, N, setup, root)
%MAKE_FIGURE  Run one panel of Fig. 3 or Fig. 4 of the paper, plot it, save curves.
%   stats = make_figure(name, beta, target_db, curves, dur_s, N, setup, root)
%
%   name       e.g. 'fig3a': writes figures/fig3a.png and results/fig3a.csv
%   beta       noise shape the filters assume (2 or 1)
%   target_db  target misalignment of the panel; selects the knobs of Tables 1-2
%   curves     C x 2 cell array of {variant, I}
%   dur_s      length of the run in seconds
%   N          number of realizations (the paper uses 100)
%   setup      struct with h, a, v_u, snr_db, beta_star, fs, seed
%
%   All curves of a panel see the same input and noise realizations.

T = round(dur_s * setup.fs);
rng(setup.seed);
x = af.ar1_input(T, N, setup.a, setup.v_u);
v_eta = af.noise_var_for_snr(setup.h, setup.snr_db, setup.a, setup.v_u);
eta = af.gg_rnd([T N], sqrt(v_eta), setup.beta_star);

C = size(curves, 1);
m_db = zeros(T, C);
labels = cell(1, C);
stats = struct('panel', {}, 'variant', {}, 'beta', {}, 'n_iter', {}, 'target_db', {}, ...
               'floor_db', {}, 't_conv', {}, 'run_s', {});
for c = 1:C
  [variant, n_iter] = curves{c, :};
  k = af.paper_params(beta, target_db, n_iter);
  opt = k.(variant);
  opt.sigma_eta = sqrt(v_eta);                  % the filters know v_eta, Sec. 4
  opt.n_iter = n_iter;
  tic;
  out = af.run_filter(setup.h, x, eta, variant, beta, opt);
  run_s = toc;
  t_sec = out.t / setup.fs;
  m_db(:, c) = 10 * log10(out.mis);
  [floor_db, t_conv] = curve_stats(out.mis, t_sec);
  if strcmp(variant, 'SG')
    labels{c} = 'SG';
  else
    labels{c} = sprintf('%s, I = %d', variant, n_iter);
  end
  stats(c) = struct('panel', name, 'variant', variant, 'beta', beta, 'n_iter', n_iter, ...
                    'target_db', target_db, 'floor_db', floor_db, 't_conv', t_conv, ...
                    'run_s', run_s);
  fprintf('  %-6s %-10s floor %6.2f dB (target %d)   within 1 dB of floor at %6.3f s   [%5.1f s]\n', ...
          name, labels{c}, floor_db, target_db, t_conv, run_s);
end

% --- plot, in the style of the paper's Figs. 3 and 4 --------------------------
colors = struct('SG', [0 0 0.9], 'fKF', [0.85 0 0], 'sKF', [0.7 0 0.7], ...
                'vKF', [0.4 0.4 0.4], 'KF', [0 0.65 0.75]);
fig = figure('visible', 'off');
hold on;
for c = 1:C
  if curves{c, 2} == 0
    style = '-';
  else
    style = '--';
  end
  plot(t_sec, m_db(:, c), style, 'Color', colors.(curves{c, 1}), 'LineWidth', 1);
end
plot([0 dur_s], [target_db target_db], ':k');
hold off;
grid on;
box on;
xlim([0 dur_s]);
ylim([target_db - 5, 2]);
xlabel('t [sec]');
ylabel('mean misalignment [dB]');
title(sprintf('\\beta = %d, \\beta^* = %.1f, target %d dB, N = %d', ...
              beta, setup.beta_star, target_db, N));
legend([labels {'target'}], 'Location', 'northeast');
print(fig, fullfile(root, 'figures', [name '.png']), '-dpng', '-r150');
close(fig);

% --- curves, thinned to about 2000 rows --------------------------------------
step = max(1, floor(T / 2000));
idx = (1:step:T)';
fid = fopen(fullfile(root, 'results', [name '.csv']), 'w');
fprintf(fid, '%% %s: mean misalignment in dB, beta = %d, target %d dB, N = %d, every %d samples\n', ...
        name, beta, target_db, N, step);
fprintf(fid, '%% t_sec');
for c = 1:C
  fprintf(fid, ', %s_I%d', curves{c, 1}, curves{c, 2});
end
fprintf(fid, '\n');
fprintf(fid, [repmat('%.6g, ', 1, C) '%.6g\n'], [t_sec(idx) m_db(idx, :)]');
fclose(fid);
end
