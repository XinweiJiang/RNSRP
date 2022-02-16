function [eigvector_new] = getPresentation_reconstruct(W,W1,train_x,option,temp_train_x)
 ratio=option.ratio;
[retrain_x]=data_RE(W,train_x,temp_train_x);
[retrain_x1]=data_RE(W1,train_x,temp_train_x);
        
[D,N] = size(train_x);
%---------------get SL-------------------------%
SL = train_x*(eye(N)-W-W'+W*W')*train_x';
SL1 = retrain_x*(eye(N)-W-W'+W*W')*retrain_x';
SL2 = train_x*(eye(N)-W1-W1'+W1*W1')*train_x';
SL3 = retrain_x1*(eye(N)-W1-W1'+W1*W1')*retrain_x1';
%---------------get ST-------------------------%
meanx   = mean(train_x,2);
ST  = zeros(D,D);
for i = 1:N
     ST = ST+ (train_x (:,i)-meanx)*(train_x(:,i)-meanx)';
end
ST=ST;
SL=(1-ratio)*(SL+SL1)+ratio*(SL2+SL3);

% SL(find(isnan(SL)==1)) = 0;
% ST(find(isnan(ST)==1)) = 0;
[eigvector, eigvalue] = eig(ST,SL);
eigvalue = diag(eigvalue);
[junk, index] = sort(-eigvalue);
eigvalue = eigvalue(index);
eigvector_new = eigvector(:,index);
for tmp = 1:size(eigvector_new,2) 
     eigvector_new(:,tmp) = eigvector_new(:,tmp)./norm(eigvector_new(:,tmp));
end
end


