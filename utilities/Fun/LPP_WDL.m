function [W,D,L]=LPP_WDL(trainX,Wk,train_total)
%根据ORL数据库的训练数据矩阵trainX，计算并返回W,D,L矩阵
%Wk是求W时需确定的参数;trainX2和trainX3分别是二维和三维的训练样本矩阵


%定义W计算式中的t
% y=pdist(trainX2);
% [t,ti]=max(y);
num=train_total;
%K=Wk时，计算出的W
W=zeros(num);
W=double(W);

%计算训练集中图像之间的距离矩阵A
for j=1:num
    for k=j:num
        z=(trainX(j,:)-trainX(k,:))*(trainX(j,:)-trainX(k,:))';
        A(j,k)=z;
        A(k,j)=z;
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 计算locality preserving structure的矩阵W、D、L
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%对A的每一行中的值进行排序，找出每一行的第(Wk+1)个最小的值，存入一维数组E中
for j=1:num
    for k=1:num
        B(k)=A(j,k);
    end
    for a=1:num-1
        for b=(a+1):num
            if (B(a)>B(b))
                y=B(a);
                B(a)=B(b);
                B(b)=y;
            end
        end
    end
    E(j)=B(Wk+1);
    t(j)=sum(B(1:Wk+1))/Wk^2;
    if (t(j)==0)
        t(j)=sum(B(1:Wk+2))/Wk^2;
    end
end
    
%根据得到的E和距离矩阵A，给矩阵W赋值
for j=1:num
    for k=1:num
        if (k~=j)
            if (A(j,k)<=E(j))
                 W(j,k)=exp(-A(j,k)/t(j));
                 W(k,j)=W(j,k);
            end
        end
    end
end
   
%计算D
D=zeros(num);
D=double(D);
for i=1:num
    D(i,i)=sum(W(i,:));
end

%计算L
L=zeros(num);
L=double(L);
L=D-W;