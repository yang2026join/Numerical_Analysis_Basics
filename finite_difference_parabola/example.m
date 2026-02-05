
%%

% Burden book, Chap 12.3, Ex 1

% Specify the PDE, BCs, ICs
a = 2;
L = 1;
T = 1;
f = @(x) sin(pi*x);
g = @(x) 0;
fdd = @(x) -pi^2 * sin(pi*x);

% Set the grid
M = 20;
N = 50;

% Apply finite-difference method
[w, x, t] = solver(a, f, g, fdd, M, N, L, T);

% Visualize
figure;
imagesc(x, t, w);
xlabel('x');  ylabel('t')
axis xy