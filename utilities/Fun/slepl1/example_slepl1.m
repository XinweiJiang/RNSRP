clear, clc;
% This is an example for applying SLEPl1 to solve
%           l1 ball constrained smooth convex optimization
% 
%               min         g(x)
%               s.t.        ||x||_1 <= z
%
%  In this example, we assume that g(x)=1/2 ||Ax-b||^2
%  By redefining the handles (line 35-38)
%            g-            g(x)
%            gprime-       g'(x)
%  it can handle general l1 ball constrained smooth convex optimization
load temp_Ay;
root=pwd;
addpath([root '/slepl1']);

m=90;  n=89;

% for reproducibility
randNum=1;

% ---------------------- generate random data ----------------------
randn('state',(randNum-1)*3+1);
%A=randn(m,n);

randn('state',(randNum-1)*3+2);
xOrin=randn(n,1);

randn('state',(randNum-1)*3+3);
noise=randn(m,1);
%y=A*xOrin + noise*0.01;
% ---------------------- generate random data ----------------------

%------------------ set the handles for g(.) and g'(.) -------------
g=@(x) 0.5* (A*x-y)' * (A*x-y);
gprime=@(x) A'*(A*x -y);
%------------------ set the handles for g(.) and g'(.) -------------

%------------------ set value for z and initialize x0 --------------
z=100;
x0=-gprime(zeros(n,1));  x0=z/sum(abs(x0)) * x0;
%------------------ set value for z and initialize x0 --------------

%------------------ run slep1 for solving 
%                          min g(x) 
%                          s.t. ||x||_1 <=z
[x, status]=slepl1(n, g, gprime, z, x0, 1, 1e-3);
plot(x)