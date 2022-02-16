function [x, status]=slepl1(n, g, gprime, z, x0, gamma, xtol, maxIter)
%
% This function solves the 
%           l1 ball constrained smooth convex optimization
%
%               min         g(x)
%               s.t.        ||x||_1 <= z
%
% Inputs: (values in * * are default values)
%  n-         the dimension of x
%  g-         a function handle that has input x and return g(x)
%  gprime-    a function handle that has input x and return g'(x)
%  z -        upper bound of ||x||_1
%  x0-        an initial guess of the solution
%             *x0=-z* gprime(zeros(n,1)) /sum(abs( gprime(zeros(n,1)) ) );*
%  gamma-     an initial guess of L_g, the Lipschitz gradient of g(.)
%             *gamma=1*
%  xtol-      the relative gap betwwen adjacent solutions
%             *tol=1e-3*
%  maxIter-   maximum number of iterations
%             *maxIter=10000*
%
% Outputs:
%  x-         the obtained solution
%  status-    1: xtol is satisfied
%             0: xtol is not satisfied 
%
% Note: 
%       1) please make sure that g(.) is smooth convex.
%       2) slepl1 does not check the correctness of g(.) and g'(.)
%       3) for Lasso, slepl1lasso is much efficient
%
% For more information on SLEPl1, 
% please refer to our paper:
% Jun Liu, Jieping Ye, and Rong Jin
% Sparse Learning with Euclidean Projection onto l_1 Ball
%
% Written by Jun Liu on Sep. 28, 2008
% Please contact: jliu86@asu.edu for any problem

if (nargin <4)
    error('\n Inputs: n, g, gprime and z should be specified! \n');
end

if (z<0)
    error('\n z should be nonnegative!\n');
end

% set default values if not specified
if (nargin<8)
    maxIter=10000;
    if (nargin<7)
        xtol=1e-3;
        if (nargin<6)
            gamma=1;
            if (nargin<5)
                x0=-gprime(zeros(n,1));   x0_abs=sum(abs(x0));
                if (x0_abs==0)
                    x=zeros(n,1); status=0; %x0 is the solution
                    return;
                else
                    x0=z/x0_abs * x0;
                end
            end
        end
    end
end

% assign xp and x with x0
x=x0; xp=x0; xxp=zeros(n,1);
tp=0; t=1;
lambda0=0; status=0; 

for iterStep=1:maxIter
    % --------------------------- step 1 ---------------------------
    % compute search point s based on xp and x (with alpha)
    alpha=(tp-1)/t;    s=x + alpha* xxp;
    
    % --------------------------- step 2 ---------------------------
    % line search for gamma and compute the new approximate solution x
    
    % compute the gradient and function value at s
    gprime_s=gprime(s);     g_s=g(s);
    
    % store x to xp
    xp=x;
    
    while (1)
        % let s walk in a step in the antigradient of s to get v
        % and project v onto the l1 ball
        v=s-gprime_s/gamma;
        [x, lambda, zf_step]=eplb(v, n, z, lambda0);
        lambda0=lambda;
        
        % compute the function value at the new approximate solution x
        g_x=g(x);
                
        v=x-s;
        % r_sum= 0.5* ||x-s||^2
        r_sum=0.5* v'*v; 
        % l_sum= g(x) -g(s) - <g'(s), x-s>
        l_sum=g_x - g_s - gprime_s'* v;
        
        % the condition is l_sum <= gamma * r_sum
        if(l_sum <= r_sum * gamma)
            break;
        else            
            gamma=max(2*gamma, l_sum/r_sum);
        end
    end
    
    % --------------------------- step 3 ---------------------------
    % update t and tp, and check whether converge
    tp=t; t= (1+ sqrt(4*t*t +1))/2; 
    
    xxp=x-xp;    norm_xp=sqrt(xp'*xp);    norm_xxp=sqrt(xxp'*xxp);
    if ( norm_xxp < max(norm_xp,1) * xtol)
        status=1; 
        break;
    end
end