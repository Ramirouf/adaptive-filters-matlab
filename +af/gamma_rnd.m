function g = gamma_rnd(a, sz)
%GAMMA_RND  Gamma(a, 1) random numbers, without the Statistics Toolbox.
%   g = af.gamma_rnd(a, sz) returns an array of size sz.
%
%   Marsaglia & Tsang, "A simple method for generating gamma variables", ACM TOMS
%   26(3), 2000. For a < 1 it uses Gamma(a) = Gamma(a + 1) * U^(1/a).
%   (randg and gamrnd would do this, but in MATLAB they need a toolbox.)

if a < 1
  g = af.gamma_rnd(a + 1, sz) .* rand(sz) .^ (1 / a);
  return;
end

d = a - 1 / 3;
c = 1 / sqrt(9 * d);
g = zeros(sz);
todo = (1:prod(sz))';
while ~isempty(todo)
  n = numel(todo);
  z = randn(n, 1);
  v = (1 + c * z) .^ 3;
  u = rand(n, 1);
  ok = v > 0 & log(u) < 0.5 * z .^ 2 + d - d * v + d * log(max(v, realmin));
  g(todo(ok)) = d * v(ok);
  todo = todo(~ok);
end
end
