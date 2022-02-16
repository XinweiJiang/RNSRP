clear, clc;

% This is an example for solving Lasso by SLEPl1 
% 
%  min     1/2 ||A x - b||^2
%  s.t.    ||x||_1 <=z
%
%  A is a unstructured dense matrix

root=pwd;
addpath([root '/slepl1']);

m=1000;  n=10000;

% for reproducibility
randNum=1;

% ---------------------- generate random data ----------------------
randn('state',(randNum-1)*3+1);
A=randn(m,n);

randn('state',(randNum-1)*3+2);
xOrin=randn(n,1);

randn('state',(randNum-1)*3+3);
noise=randn(m,1);
y=A*xOrin + noise*0.01;
% ---------------------- generate random data ----------------------

% ----------------------------- SLEPl1 -----------------------------
z=1000;
x0=A'*y; x0=z/sum(abs(x0)) * x0;
tic;
[x, reg_rho, fun_val, iterStep, status]=slepl1_lasso(A, m, n, y, z, x0, 1, 1e-3);
t=toc;
% ----------------------------- SLEPl1 -----------------------------

figure;
plot(reg_rho)
title('A: dense')
xlabel('iteration (i)');
ylabel('regularization: \lambda_i\gamma_i');

figure;
plot(fun_val)
title('A: dense')
xlabel('iteration (i)');
ylabel('objective function value');

figure;
plot(iterStep)
title('A: dense')
xlabel('iteration (i)');
ylabel('number of inner iterations for finding the root');