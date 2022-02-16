function [eigvector_new] = getPresentation(W,train_x)
[D,N] = size(train_x);
%---------------get SL-------------------------%
SL = train_x*(eye(N)-W-W'+W*W')*train_x';
%---------------get ST-------------------------%
meanx  = mean(train_x,2);
ST  = zeros(D,D);
for i = 1:N
    ST = ST+ (train_x(:,i)-meanx)*(train_x(:,i)-meanx)';
end
% SL=max(SL,SL');
% ST=max(ST,ST');
[eigvector, eigvalue] = eig(ST,SL);
eigvalue = diag(eigvalue);
[junk, index] = sort(-eigvalue);
eigvalue = eigvalue(index);
eigvector_new = eigvector(:,index);
for tmp = 1:size(eigvector_new,2) 
     eigvector_new(:,tmp) = eigvector_new(:,tmp)./norm(eigvector_new(:,tmp));
end

end

