
%% Specify the boundary-value problem
% adapted from Burden book, Chap. 12.1, Ex. 1

f = @(x,y) 0*x+0*y;
a = 0;
b = 0.5;
c = 0;
d = 0.5;
bc_ay = @(y) 25*y;
bc_by = @(y) 25 + (400-25)/0.5*y;
bc_xc = @(x) 25*x;
bc_xd = @(x) 25 + (400-25)/0.5*x;
N = 6;
M = 6;
maxit = 100;
TOL = 10e-5;

%% Visualize the linear system

[A, b, b1, b2, b3] = finite_difference_elliptic_linear_system(f, a, b, c, d, bc_ay, bc_by, bc_xc, bc_xd, N, M);

figure;
imagesc(A);  colorbar;  title('Matrix A')
axis square
figure;
subplot(4,1,1);  imagesc(b');  colorbar;  title('b');  set(gca, 'YTick', []);
subplot(4,1,2);  imagesc(b1');  colorbar;  title('b_1');  set(gca, 'YTick', []);
subplot(4,1,3);  imagesc(b2');  colorbar;  title('b_2');  set(gca, 'YTick', []);
subplot(4,1,4);  imagesc(b3');  colorbar;  title('b_3');  set(gca, 'YTick', []);
%set(gcf, 'Position', [100 100 600 60]);

%% Apply finite-difference method and visualize the solution

[w, x, y] = finite_difference_elliptic(f, a, b, c, d, bc_ay, bc_by, bc_xc, bc_xd, N, M, maxit, TOL);

figure;
imagesc(x, y, w);  colorbar;
xlabel('x');  ylabel('y');
title('Steady-state Temperature Distribution');
axis xy
axis square
