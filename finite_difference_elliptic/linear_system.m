function [A, b, b1, b2, b3] = linear_system(f, a, b, c, d, bc_ay, bc_by, bc_xc, bc_xd, N, M)

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

b1 = f(x(2:N), flip(y(2:M)'));
b1 = reshape(b1.', 1, [])'*(-h^2);

b2 = [w_iM zeros(1,(N-1)*(M-3)) w_i0]';
b2 = lambda*b2;

b3 = [flip(w_0j(2:M)) zeros(M-1,N-3) flip(w_Nj(2:M))];
b3 = reshape(b3.', 1, [])';

b = b1 + b2 + b3;

end