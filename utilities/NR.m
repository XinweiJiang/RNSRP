function [Z,C] = DNRC_41(X, Y, param)
mu = param.beta;
 lambda = param.lamuda;
[~,n] = size(X);
m = size(Y,2);
tol = 1e-5;
maxIter = 5;
I=eye(n);
% initialization
Z = zeros(n,m);
C = zeros(n,m);
delta = zeros(n,m);

XTX = X'*X;
XTY = X'*Y;
iter = 0;

% pre-computation
temp_X = pinv((1+ lambda)*XTX+mu/2*eye(n));%对于奇异矩阵或者非方阵，并不存在逆矩阵，但可以使用pinv(A)求其伪逆

while iter<maxIter
    iter = iter + 1;
    
    Zk = Z;
    Ck = C;
    
    % update c
    C = temp_X*(XTY+mu/2*Z+delta/2);
    
    % update z
    z_temp = C-delta/mu;
    Z = max(0,z_temp);
    
    leq1 = Z-C;
    leq2 = Z-Zk;
    leq3 = C-Ck;
    stopC1 = max(norm(leq1,'fro'),norm(leq2,'fro'));
    stopC = max(stopC1,norm(leq3,'fro'));
    %     disp(stopC)
    
    if stopC<tol || iter>=maxIter
        break;
    else
        % update delta
        delta = delta + mu*leq1;
    end
end