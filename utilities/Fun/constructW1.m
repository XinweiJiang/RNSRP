%===Graph Construction based on Sparse Representation(SR)=================%
%
%Input:
%    X                - data, each colummn denotes a simple.
%    options.mode     - 'standard' or 'extend' mode. default='extend'.
%    options.lambda   - regularized parameter for extended mode, default =1.
%Output:
%    Affinity weight matrix W.
%
%===Version 1.0, Written by Lishan Qiao, 2009.=============================

function [W] = constructW1(X,options)

path(path, './Fun/l1magic/Optimization');
[D,N] = size(X); W = zeros(N,N);

if strcmp(options.mode,'standard')
    %%---calculate sparse reconstruction weights---%%
    for i=1:N
        y=[X(:,i);1];%% for invariance to translations
        A=[X(:,setdiff(1:N,i)); ones(1,N-1)];
        w0 = A\y;
        wi = l1eq_pd(w0, A, [], y, 1e-3);%%%%%%%%%---!!!!!!!!!!---%%%%%%%%%
        W(setdiff(1:N,i),i) = wi;
        fprintf('No.%d sample is done\n',i)
    end
elseif strcmp(options.mode,'extend')
    if ~isfield(options,'lambda')
        options.lambda=1;
    end
    for i=1:N
        y=[X(:,i);1];%% for invariance to translations
        A=[[X(:,setdiff(1:N,i)); ones(1,N-1)] [eye(D)/options.lambda; zeros(1,D)]];
        w0 = A\y;
        wi = l1eq_pd(w0, A, [], y, 1e-3);%%%%%%%%%---!!!!!!!!!!---%%%%%%%%%
        W(setdiff(1:N,i),i) = wi(1:N-1);
        fprintf('No.%d sample is done\n',i)
    end
end