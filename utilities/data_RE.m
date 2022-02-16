function [re_x,W_1]=data_RE(W,train_x,temp_train_x)
[D,N]=size(train_x);
re_x = zeros(D,N);
W_1=zeros(N,N);
[tmp, ind] = sort(W,'descend');
for i = 1:N 
     w=zeros(N,1);
     nf=size(find(W(:,i)>0),1);
     w(ind(1:nf,i))=W(ind(1:nf,i),i);
     W_1(ind(1:nf,i),i)=W(ind(1:nf,i),i);
     temp = w./sum(w);
     re_x(:,i) = temp_train_x*temp;
end