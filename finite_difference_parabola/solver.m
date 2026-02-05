function [w, x, t] = solver(a, f, g, fdd, M, N, L, T)
%Apply finite-difference method to approximate the solution to the linear
% second-order boundary-value problem: u_tt = a^2 * u_xx
% with boundary conditions u(0, t) = u(1, t) = 0 for t >= 0
% and initial conditions u(x, 0) = f(x), u_x(x, 0) = g(x), for 0 <= x <= L
%--------------------------------------------------------------------------
% Input
% a: parameter in the PDE
% f, g: function handle of f(x), g(x) in the PDE; fdd: f''(x)
% l, T: range of x, maximum time
% M, N: number of grid lines for x-axis, t-axis
%--------------------------------------------------------------------------
% Output
% w: approximations w_i,j to u(x_i, t_j)
% for i = 0, 1, 2, ..., M and j = 0, 1, 2, ..., N
%--------------------------------------------------------------------------

% time step
h = L/M;
k = T/N;
lambda = a*k/h;
x = linspace(h, L-h, M-1);    % interior points (i = 1, 2, ..., M-1)
t = linspace(k, T, N);    % time points

% form the transitional matrix
D = sparse(1:M-1, 1:M-1, 2*(1-lambda^2)*ones(1,M-1));
L = sparse(2:M-1, 1:M-2, lambda^2*ones(1,M-2), M-1, M-1);
U = sparse(1:M-2, 2:M-1, lambda^2*ones(1,M-2), M-1, M-1);
A = D + L + U;

% iterate over time
w = zeros(M-1, N);
w(:, 1) = f(x);    % j = 0
w(:, 2) = w(:, 1) + k*g(x)' + a^2*k^2*fdd(x)'/2;    % j = 1
for j = 2:(N-1)
    w(:, j+1) = A*w(:, j) - w(:, j-1);
end

end