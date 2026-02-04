function [w, x, y] = solver(f, a, b, c, d, bc_ay, bc_by, bc_xc, bc_xd, N, M, maxit, TOL)
%FINITE_DIFFERENCE_ELLIPTIC
% to solve the boundary problem of elliptic equation:
% Δu = f(x, y) on rectangle area R = {(x, y) | a < x < b, c < y < d}
% with boundary condition u(x, y) = bc(x, y) for (x, y) ∈ Σ(R),
%--------------------------------------------------------------------------
% Input
% f: function handle for f(x, y)
% bc_ay, bc_by, bc_xc, bc_xd: function handles for boundary conditions
% a, b, c, d: range of the rectangle area
% N, M: create an (N + 1) × (M + 1) mesh grid
% maxit: maximum number of iterations
% TOL: tolerance for iterative method
%--------------------------------------------------------------------------
% Output
% x, y: mesh points, a + ih and b + jk, for i = 0, 1, ..., N + 1, and j = 0, 1, ..., M + 1
% w: approximated values at the mesh points (x_i, y_j)
%--------------------------------------------------------------------------

%% mesh points

x = linspace(a, b, N+1);    % x_0, x_1, ..., x_N-1, x_N
y = linspace(c, d, M+1);    % y_0, y_1, ..., y_M-1, y_M
h = (b - a) / N;
k = (d - c) / M;
lambda = h^2/k^2;
mu = 2 * (lambda + 1);

%% boundary values

% Note: avoid repeated assignment of values at the corners
w_0j = bc_ay(y)';    % j = 0, 1, ..., M
w_Nj = bc_by(y)';    % j = 0, 1, ..., M
w_i0 = bc_xc(x(2:end-1));    % i = 1, 2, ..., N-1
w_iM = bc_xd(x(2:end-1));    % i = 1, 2, ..., N-1

%% form matrix A

disp('Form matrix A.')

i = 1:N-1;  j = 1:N-1;  v = -lambda*ones(1,N-1);
Lambda = sparse(i, j, v);

i1 = 1:N-1;  j1 = 1:N-1;  v1 = mu*ones(1,N-1);
i2 = 1:N-2;  j2 = 2:N-1;  v2 = -ones(1,N-2);
i3 = 2:N-1;  j3 = 1:N-2;  v3 = -ones(1,N-2);
Mu = sparse([i1 i2 i3], [j1 j2 j3], [v1 v2 v3]);

A = [Mu Lambda];
for i = 2:(M-2)
    A = [A zeros((i-1)*(N-1), N-1); zeros(N-1, (i-2)*(N-1)) Lambda Mu Lambda];
end
A = [A; zeros(N-1, (M-3)*(N-1)) Lambda Mu];

%% form vector b

disp('Form vector b.')

b1 = f(x(2:N), flip(y(2:M)'));
b1 = reshape(b1.', 1, [])'*(-h^2);

b2 = [w_iM zeros(1,(N-1)*(M-3)) w_i0]';
b2 = lambda*b2;

b3 = [flip(w_0j(2:M)) zeros(M-1,N-3) flip(w_Nj(2:M))];
b3 = reshape(b3.', 1, [])';

b = b1 + b2 + b3;

%% apply iterative method

disp('Apply iterative method to solve the linear system.')

% apply Jacobian iterative technique with SOR
weight = 4 / (2 + sqrt(4 - (cos(pi/M) + cos(pi/N))^2));    % Burden book p740
x0 = randn((M-1)*(N-1), 1);
w = SOR(A, b, x0, TOL, maxit, weight);
%size(A), size(b)
% For small systems, can solve directly.
%w = A\b;

% organize as 2-d mesh grid
w = flip(reshape(w', M-1, [])', 1);

% add the boundary values
w = [w_i0; w; w_iM];
w = [w_0j w w_Nj];

end

%% Jacobian with SOR iterative technique
function y = SOR(A, b, x0, TOL, maxit, w)
%--------------------------------------------------------------------------
% Input
% n:     number of equations and unknowns
% A,b:   in Ax = b
% x0:    initial guess
% TOL:   tolerance
% maxit: maximum number of iterations
% w:     weight
%--------------------------------------------------------------------------
% Output
% approximation of x in each iteration
%--------------------------------------------------------------------------

n = numel(b);    % number of equations and unknowns
x1 = zeros(n, 1);

% loop over iterations
for k = 1:maxit

    % loop over each component of x
    for i = 1:n
        x1(i) = (1 - w) * x0(i) + w / A(i, i) * (b(i)...
            - dot(A(i, 1:i-1), x1(1:i-1))...
            - dot(A(i, i+1:n), x0(i+1:n)) );
    end

    % criterion of termination
    res = norm(x1 - x0);
    fprintf("Iteration %d. Residual = %.4f\n", k, res);
    if res < TOL
        break;
    else
        x0 = x1;
    end

end

y = x1;
end