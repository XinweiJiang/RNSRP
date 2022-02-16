function [W,D,L]=LPP_WDL(trainX,Wk,train_total)
%����ORL���ݿ��ѵ�����ݾ���trainX�����㲢����W,D,L����
%Wk����Wʱ��ȷ���Ĳ���;trainX2��trainX3�ֱ��Ƕ�ά����ά��ѵ����������


%����W����ʽ�е�t
% y=pdist(trainX2);
% [t,ti]=max(y);
num=train_total;
%K=Wkʱ���������W
W=zeros(num);
W=double(W);

%����ѵ������ͼ��֮��ľ������A
for j=1:num
    for k=j:num
        z=(trainX(j,:)-trainX(k,:))*(trainX(j,:)-trainX(k,:))';
        A(j,k)=z;
        A(k,j)=z;
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% ����locality preserving structure�ľ���W��D��L
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%��A��ÿһ���е�ֵ���������ҳ�ÿһ�еĵ�(Wk+1)����С��ֵ������һά����E��
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
    
%���ݵõ���E�;������A��������W��ֵ
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
   
%����D
D=zeros(num);
D=double(D);
for i=1:num
    D(i,i)=sum(W(i,:));
end

%����L
L=zeros(num);
L=double(L);
L=D-W;