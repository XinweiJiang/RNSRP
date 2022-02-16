clear, clc;

% This is an example for solving pathwise Lasso by SLEPl1 
%       with the function: 
%
% 
%  min     1/2 ||A x - b||^2
%  s.t.    ||x||_1 <=z
%
%  A is a unstructured dense matrix
%  the pathwise solutions corresponding to a series of z in zeta are
%  computed with the "warm" start technique

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

zeta=[0.1, 10, 100, 300, 350, 400, 450, 500];
[X, status]=slepl1pathlasso(A, m, n, y, zeta);