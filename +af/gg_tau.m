function tau = gg_tau(sigma, beta)
%GG_TAU  Noise constant tau of the robust filters, Eq. (69).
%   tau = af.gg_tau(sigma, beta) = (sigma * kappa(beta))^beta / beta
%
%   It enters the gain multiplier as alpha_t = 1 / (tau |e_t|^(2-beta) + s_t), Eq. (75).
%   beta = 2: tau = sigma^2 = v_eta, so alpha_t = 1 / (v_eta + s_t), the Kalman filter.
%   beta = 1: tau = sigma / sqrt(2). Eq. (80) prints sigma * sqrt(2), 2x larger: that
%   is what (69) gives if kappa(1) were sqrt(2), but Eq. (61) gives 1/sqrt(2). The
%   general form (69) is used here.

tau = (sigma .* af.gg_kappa(beta)) .^ beta ./ beta;
end
