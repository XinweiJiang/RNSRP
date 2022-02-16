function [ x_training,xClean_training,y_training,x_testing,y_testing] = ChooseRSdata_zxx(x,y ,nTrEachClass, nSeed,radius)
% Randomly split HSI data to training and testing data 
% Input:
% dataSetName: HSI data downloaded from http://alweb.ehu.es/ccwintco/index.php?title=Hyperspectral_Remote_Sensing_Scenes
% nTrEachClass: 
%         >1: the number of training data chosen from each class
%         0-1:  the percentage of data to be training data chosen from each class
% nSeed: Random function seed
% 
% Output:
% x_training: NxD training data where each row is a sample
% y_training: Nx1 labels correponding to the training data
% x_test: NNxD testing data where each row is a sample
% x_test: NNx1 labels correponding to the testing data
addpath(genpath('utilities\.'));
xClean = x;
if radius > 0
    for i=1:size(x,3)
        x(:,:,i) = imfilter(x(:,:,i),fspecial('average',radius));
    end
end


Dim = size(x);
ind = tabulate(y(:));
nClass = size(ind,1)-1;    %except label=0 
xClean=reshape(xClean,Dim(1)*Dim(2),Dim(3));
x=reshape(x,Dim(1)*Dim(2),Dim(3));
y=reshape(y,Dim(1)*Dim(2),1);
 
xClean = sgpNormalize( xClean);
x = sgpNormalize( x);

 %
y_training=[];%
x_training=[];%
y_testing=[];%
x_testing=[];%
index_training=[];
index_testing=[];
nMaxPercent = 0.6;      %max percentage of data chosen from each class
nTrEachInit = nTrEachClass;
 for i=1:nClass
     nTrEachClass = nTrEachInit;
     % choose training data according to the percentage from each class
     if nTrEachClass < 1
         if nTrEachClass <= 0
             error('Please check input: nTrEachClass');
         end
          A= find(y==ind(i+1,1));
        nTrEachClass = round(length(A)*nTrEachClass);
        %min number of data chosen from each class
        %which void that small number of testing data cause problem for KNN
        %in graph embedding
       if nTrEachClass < 5
           nTrEachClass = 5;
       end

       rng((i+nSeed)*10,'twister');
       randomorder = randperm(length(A));
       index_training = [index_training;A(randomorder(1:nTrEachClass))];
       index_testing = [index_testing;A(randomorder(nTrEachClass+1:end))];%
     else        
         % choose training data according to the number from each class
        A= find(y==ind(i+1,1));
        %max percentage of data chosen from each class
        %which void that all data from some class to be training data
        %leading to no testing data
        if round(length(A)*nMaxPercent) < nTrEachClass    
            nTrEachClass = round(length(A)*nMaxPercent);
        end
       rng((i+nSeed)*10,'twister');
       randomorder = randperm(length(A));
       index_training = [index_training;A(randomorder(1:nTrEachClass))];
       index_testing = [index_testing;A(randomorder(nTrEachClass+1:end))];
     end
 end

x_training = [x(index_training,:)];
y_training = y(index_training,:);
x_testing = [x(index_testing,:)];
y_testing = y(index_testing,:);
xClean_training =[xClean(index_training,:)]; 

 clear x y x_clean;
 
 % 
if nClass == 2
    y_training(y_training==1) = -1; y_training(y_training==2) = 1;
    y_test(y_test==1) = -1; y_test(y_test==2) = 1;
else
    y_training = smgpTransformLabel( y_training );
    y_testing = smgpTransformLabel( y_testing );
    y_training = smgpTransformLabel( y_training );
    y_testing = smgpTransformLabel( y_testing );
end


end