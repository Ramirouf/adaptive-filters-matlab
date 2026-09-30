function e = gg_rnd(sz, sigma, beta)
%GG_RND  Zero-mean generalized Gaussian noise, Eqs. (59)-(62).
%   e = af.gg_rnd(sz, sigma, beta) returns an array of size sz with standard
%   deviation sigma and shape beta (2 = Gaussian, 1 = Laplace, < 1 heavy-tailed).
%
%   If G ~ Gamma(1/beta, 1), then |e| = c G^(1/beta), with c = sigma kappa(beta),
%   has the right magnitude distribution; the sign is drawn independently.

c = sigma * af.gg_kappa(beta);
mag = c * af.gamma_rnd(1 / beta, sz) .^ (1 / beta);
e = mag .* (2 * (rand(sz) < 0.5) - 1);
end
