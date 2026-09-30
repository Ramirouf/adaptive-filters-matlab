function x = ar1_input(T, N, a, v_u)
%AR1_INPUT  AR(1) input of Sec. 4: x_t = -a x_{t-1} + u_t, u_t ~ N(0, v_u).
%   x = af.ar1_input(T, N, a, v_u) returns T x N, one realization per column.
%   Each column starts from the stationary distribution, variance v_u / (1 - a^2).

var_x = v_u / (1 - a ^ 2);
x0 = sqrt(var_x) * randn(1, N);
u = sqrt(v_u) * randn(T, N);
x = filter(1, [1 a], u, -a * x0);
end
