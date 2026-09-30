function v_eta = noise_var_for_snr(h, snr_db, a, v_u)
%NOISE_VAR_FOR_SNR  Noise variance v*_eta that gives the output SNR of Eq. (89).
%   v_eta = af.noise_var_for_snr(h, snr_db, a, v_u)
%
%   SNR = 10 log10(E[(x_t' h)^2] / v*_eta). For the AR(1) input of Sec. 4,
%   E[(x_t' h)^2] = h' R h with R(i,j) = var_x (-a)^|i-j| and var_x = v_u / (1 - a^2).

h = h(:);
M = numel(h);
var_x = v_u / (1 - a ^ 2);
R = var_x * toeplitz((-a) .^ (0:M-1));
v_eta = (h' * R * h) / 10 ^ (snr_db / 10);
end
