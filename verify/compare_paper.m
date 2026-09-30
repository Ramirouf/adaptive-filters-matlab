function ok = compare_paper(root)
%COMPARE_PAPER  Compare our Figs. 3a and 3b with the markers of the published figures.
%   ok = compare_paper(root) reads data/paper_fig3_markers.csv (read off the paper's
%   figures by data/digitize_fig3.py) and results/fig3a.csv, results/fig3b.csv (written
%   by run_all), and prints the gap at every marker. Around each marker time, our curve
%   is averaged on the linear scale over +-1/120 of the panel's length (0.1 s in Fig. 3a).
%
%   Tolerance 1 dB. The reading is good to about 0.1 dB; the rest of the gap is Monte
%   Carlo noise, ours and the paper's (each averages N runs of heavy-tailed noise).

tol = 1;
fprintf('\n== Comparison with the markers of the published Figs. 3a and 3b ==\n');

fid = fopen(fullfile(root, 'data', 'paper_fig3_markers.csv'), 'r');
c = textscan(fid, '%s %s %f %f', 'Delimiter', ',', 'CommentStyle', '%');
fclose(fid);
panel = strtrim(c{1});
variant = strtrim(c{2});
t_mk = c{3};
m_mk = c{4};

columns = {'SG', 'fKF', 'sKF', 'KF'};          % column order of results/fig3*.csv
ok = true;
worst = 0;
fprintf('  %-6s %-5s %8s %10s %10s %8s\n', 'panel', 'filter', 't [s]', 'paper [dB]', 'ours [dB]', 'gap');
for p = {'fig3a', 'fig3b'}
  name = p{1};
  file = fullfile(root, 'results', [name '.csv']);
  if ~exist(file, 'file')
    fprintf('  %s missing: run run_all first\n', file);
    ok = false;
    continue;
  end
  fid = fopen(file, 'r');
  header = fgetl(fid);
  fclose(fid);
  n_runs = str2double(regexp(header, 'N = (\d+)', 'tokens', 'once'));
  d = load(file, '-ascii');
  t = d(:, 1);
  half = t(end) / 120;
  for i = find(strcmp(panel, name))'
    col = 1 + find(strcmp(columns, variant{i}));
    sel = abs(t - t_mk(i)) <= half;
    ours = 10 * log10(mean(10 .^ (d(sel, col) / 10)));
    gap = ours - m_mk(i);
    worst = max(worst, abs(gap));
    pass = abs(gap) < tol;
    ok = ok && pass;
    fprintf('  %-6s %-5s %8.3f %10.2f %10.2f %+8.2f   %s\n', name, variant{i}, t_mk(i), ...
            m_mk(i), ours, gap, pass_str(pass));
  end
  if n_runs < 100
    fprintf('  (%s uses N = %d, not the paper''s 100: expect larger gaps)\n', name, n_runs);
  end
end

fprintf('\n  largest gap %.2f dB over %d markers (tolerance %g dB)\n', worst, numel(t_mk), tol);
if ok
  fprintf('  comparison with the paper PASSED\n');
else
  fprintf('  comparison with the paper FAILED\n');
end
end

function s = pass_str(pass)
if pass
  s = 'PASS';
else
  s = 'FAIL';
end
end
