function [d_train_x,d_test_x] = RNSRP(train_x,temp_train_x,test_x,options,maxDim)

W=creatgraphW('NRP',train_x,options);
W1=creatgraphW('SPP',train_x,options);
eigvector = getPresentation_reconstruct(W,W1,train_x,options,temp_train_x);
%  eigvector = getPresentation(W,train_x);
P =  eigvector(:,1:maxDim)';
    if ~isreal(P)
        warning('there is complex data!');
    end
 d_train_x= P*train_x;
 d_test_x = P*test_x;
end


