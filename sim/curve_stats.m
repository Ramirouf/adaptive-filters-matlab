function [floor_db, t_conv] = curve_stats(mis, t_sec)
%CURVE_STATS  Steady-state floor and convergence time of one misalignment curve.
%   [floor_db, t_conv] = curve_stats(mis, t_sec)
%
%   floor_db  mean of the last 10% of the curve (on the linear scale), in dB
%   t_conv    first time the curve, smoothed over 1% of its length, is within
%             1 dB of floor_db, in seconds

mis = mis(:);
n = numel(mis);
floor_db = 10 * log10(mean(mis(n - max(1, round(0.1 * n)) + 1:n)));
w = max(1, round(0.01 * n));
sm = conv(mis, ones(w, 1) / w, 'valid');        % sm(k) is centred on sample k + (w-1)/2
k = find(10 * log10(sm) <= floor_db + 1, 1);
t_conv = t_sec(k + floor((w - 1) / 2));
end
