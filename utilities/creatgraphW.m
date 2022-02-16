
function W = creatgraphW(type,Datatrain,options)
%
W=[]; 
[demnsion,samples]=size(Datatrain);
switch upper(type)
    
case upper('SPP')
    path(path, './Fun/l1magic/Optimization');
     W = zeros(samples,samples);
    %%---calculate sparse reconstruction weights---%%
    for i=1:samples
        [D,N] = size(Datatrain);
        ind=[];
        y=[Datatrain(:,i);1];%% for invariance to translations
        A=[[Datatrain(:,setdiff(1:samples,i)); ones(1,samples-1)] [eye(D)/10^0; zeros(1,D)]];%10^0
        w0 = A\y;
        wi = l1eq_pd(w0, A, [], y, 1e-3);
        ind=find(wi<0.00001);
        wi(ind)=0;
        W(setdiff(1:samples,i),i) = wi(1:samples-1);
        fprintf('No.%d sample is done\n',i)
    end
case upper('NRP')
        for i = 1: samples
            fprintf('NRP use:I calculate the reconstruct weight of %d,total %d datapoint\n',[i,samples]);
            temp = Datatrain;   % Da ta: d x m
            datai = temp(:, i); % Y: d x 1
            temp(:, i) = [];
            [A,~] = NR(temp, datai,options);
            B = zeros(samples, 1);
            B([1:i-1, i+1:samples]) = A;
            W(:,i) = B;
        end
    otherwise
        error('输入的创建形式有误！');
end


        
