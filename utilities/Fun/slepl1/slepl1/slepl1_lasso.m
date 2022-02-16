function [x, reg_rho, fun_val, rootStep, status]=...
    slepl1_lasso(A, m, n, y, z, x0, gamma, xtol, maxIter)
%
% slepl1_lasso: apply slepl1 to solving the lasso problem
%               with outputs such as reg_rho, fun_val, rootStep
%
%  min  1/2 || A x - y||
%  s.t. ||x||_1 <= z
%
% Inputs: (values in * * are default values)
%  A -        matrix of size m x n
%             A can be a dense matrix
%                      a sparse matrix
%                      or a DCT matrix
%  y -        response vector (of size nx1)
%  z -        upper bound of ||x||_1
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
%  x-         the obtained solution
%  rho-       the associated regularization parmameter
%  rootStep-  number of iterations in zero finding
%  status-    1: xtol is satisfied
%             0: xtol is not satisfied 
%
% For more information on SLEPl1, 
% please refer to our paper:
% Jun Liu, Jieping Ye, and Rong Jin
% Sparse Learning with Euclidean Projection onto l_1 Ball
%
% Written by Jun Liu on Sep. 28, 2008
% Please contact: jliu86@asu.edu for any problem

if (nargin <5)
    error('\n Insufficient inputs: A, y and z should be specified');
end

if (length(y) ~=m)
    error('\n Check the length of y');
end

% set default values if not specified
if (nargin<9)
    maxIter=10000;
    if (nargin<8)
        xtol=1e-3;
        if (nargin<7)
            gamma=1;
            if (nargin<6)
                x0=A'*y; x0=z/sum(abs(x0)) * x0;
            end
        end
    end
end

% assign xp and x with x0, and compute Ax
x=x0; xp=x0;
Ax=A*x; Axp=Ax; xxp=x -xp;

ATy=A'*y;  tp=0; t=1;
lambda0=0; status=0; 

for iterStep=1:maxIter
    % --------------------------- step 1 ---------------------------
    % compute search point s based on xp and x (with alpha)
    alpha=(tp-1)/t;    s=x + alpha* xxp;
    
    % --------------------------- step 2 ---------------------------
    % line search for gamma and compute the new solution x
    
    % compute the gradient (g) at s
    As=Ax + alpha* (Ax-Axp);    g=A'*As-ATy;  
    
    % copy x and Ax to xp and Axp
    xp=x;    Axp=Ax;
    
    while (1)
        % let s walk in a step in the antigradient of s to get v
        % and project v onto the l1 ball
        v=s-g/gamma;
        [x, lambda, zf_step]=eplb(v, n, z, lambda0);
        lambda0=lambda; rootStep(iterStep)=zf_step;
        
        v=x-s;  Ax=A*x;  Av=Ax -As;
        r_sum=v'*v; l_sum=Av'*Av;
        
        % the condition is ||A * v||_2^2 <= gamma * ||v||_2^2
        if(l_sum <= r_sum * gamma)
            reg_rho(iterStep)=gamma* lambda;
            break;
        else            
            gamma=max(2*gamma, l_sum/r_sum);
            %fprintf('\n gamma=%5.6f',gamma);
        end
    end
    
    % --------------------------- step 3 ---------------------------
    % update t and tp, and check whether converge
    tp=t; t= (1+ sqrt(4*t*t +1))/2; 
    
    fun_val(iterStep)=0.5* (Ax-y)'* (Ax-y);
   
    xxp=x-xp;    norm_xp=sqrt(xp'*xp);    norm_xxp=sqrt(xxp'*xxp);
    if ( norm_xxp < max(norm_xp,1) * xtol)
        status=1;
        break;
    end
end

