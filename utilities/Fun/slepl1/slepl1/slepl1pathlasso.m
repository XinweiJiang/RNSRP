function [X, status]=slepl1pathlasso(A, m, n, y, zeta, x0, gamma, xtol, maxIter)
%
% slepl1pathlasso: apply slepl1 to solving pathwise lasso
%
%  min  1/2 || A x - y||
%  s.t. ||x||_1 <= z
%
%  with z sequentially chosen from zeta
% 
%
% Inputs: (values in * * are default values)
%  A -        matrix of size m x n
%             A can be a dense matrix
%                      a sparse matrix
%                      or a DCT matrix
%  y -        response vector (of size nx1)
%  zeta -     a vector contains increasing values of z's
%  x0-        an initial guess of the solution
%             *x0=A'*y; x0=z/sum(abs(x0)) * x0;*
%  gamma-     an initial guess of the maximal eigenvalue of A'A
%             *gamma=1*
%  xtol-      the relative gap betwwen adjacent solutions
%             *tol=1e-3*
%  maxIter-   maximum number of iterations
%             *maxIter=10000*
%
% Outputs:
%  X-         the pathwise solutions corresponding to each z in zeta
%  status-    1: xtol is satisfied for all z's in zeta
%             0: xtol is not satisfied for any z in zeta
%
% For more information on SLEPl1, 
% please refer to our paper:
% Jun Liu, Jieping Ye, and Rong Jin
% Sparse Learning with Euclidean Projection onto l_1 Ball
%
% Written by Jun Liu on Sep. 28, 2008
% Please contact: jliu86@asu.edu for any problem

if (nargin <5)
    error('\n Inputs: A, m, n, y and z should be specified!\n');
end

if (length(y) ~=m)
    error('\n Check the length of y!\n');
end

k=length(zeta); status=1;
[z_value, ind]=sort(zeta);
% ensure that the values are in an increasing order

if (z_value(1) <=0)
    error('\n Values in zeta should be positive!\n');
end

% set default values if not specified
if (nargin<9)
    maxIter=10000;
    if (nargin<8)
        xtol=1e-3;
        if (nargin<7)
            gamma=1;
            if (nargin<6)
                x0=A'*y; x0_abs=sum(abs(x0));
                if (x0_abs==0)
                    X=zeros(n,k); % the solutios are zeros
                    return;
                else
                    x0=z_value(1)/x0_abs * x0;
                end
            end
        end
    end
end

for i=1:k
    z=z_value(i);

    [x, sta]=slepl1lasso(A, m, n, y, z, x0, gamma, xtol, maxIter);
    
    X(:,ind(i))=x;
    if (~sta)
        status=0;
    end

    % set x0=x, so that x0 acts as a "warm" start for the larger z
    x0=x;
end