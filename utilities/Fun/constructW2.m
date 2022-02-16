%===Graph Construction based on Sparse Representation(SR)=================%
%
%Input:
%    X                - data, each colummn denotes a simple.
%    options.mode     - 'magic' or 'slep' optimization toolbox.
%    options.epsilon  - for 'l1magic', min ||x||_1  s.t. ||y-Xw||<=epsilon, default =0.05.
%    options.k        - for 'slepl1', min ||y-Xw|| s.t. ||x||_1 <=k, default = 5.
%Output:
%    Affinity weight matrix W.
%
%===Version 1.0, Written by Lishan Qiao, 2009.=============================

function [W] = constructW2(X,options)

path(path, './Fun/l1magic/Optimization');
path(path, './Fun/slepl1/slepl1');
[D,N] = size(X); W = zeros(N,N);


if strcmp(options.mode,'magic')
    %%---calculate sparse reconstruction weights---%%
    if ~isfield(options,'epsilon')
        options.epsilon=0.05;
    end
    for i=1:N
        y=[X(:,i);1];%% for invariance to translations
        A=[X(:,setdiff(1:N,i)); ones(1,N-1)];
        w0 = A\y;
        wi = l1qc_logbarrier(w0, A, [], y, options.epsilon, 1e-3);%%%%%%%%%---!!!!!!!!!!---%%%%%%%%%
        W(setdiff(1:N,i),i) = wi;
        fprintf('No.%d sample is done\n',i)
    end
elseif strcmp(options.mode,'slep')
    if ~isfield(options,'k')
        options.k = 5;
    end
    for i=1:N
        y=[X(:,i);1];%% for invariance to translations
        A=[X(:,setdiff(1:N,i)); ones(1,N-1)];
        g=@(w) 0.5* (A*w-y)' * (A*w-y);
        gprime=@(w) A'*(A*w -y);
        z=options.k;
        w0=-gprime(zeros(N-1,1));  w0=z/sum(abs(w0)) * w0;
        [wi, status]=slepl1(N-1, g, gprime, z, w0, 1, 1e-3);%%%%%%%%%---!!!!!!!!!!---%%%%%%%%%
        W(setdiff(1:N,i),i) = wi(1:N-1);
        fprintf('No.%d sample is done\n',i)
    end
end