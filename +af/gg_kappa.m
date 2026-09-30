function k = gg_kappa(beta)
%GG_KAPPA  Shape constant of the generalized Gaussian, Eq. (61).
%   k = af.gg_kappa(beta) = sqrt(Gamma(1/beta) / Gamma(3/beta))
%
%   It converts a standard deviation into the scale of the distribution,
%   c = sigma * kappa(beta), Eq. (62). kappa(2) = sqrt(2), kappa(1) = 1/sqrt(2).

k = sqrt(gamma(1 ./ beta) ./ gamma(3 ./ beta));
end
