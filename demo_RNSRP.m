%The code is for the following work, if you find it is useful, please cite our work:
%
%[1] L. Xiong, X. Jiang, P. Lu, Z. Chen, Y. Zhang, X. Liu and Z. Cai, "Unsupervised Dimensionality Reduction for Hyperspectral Imagery via Non-negative Representation Projection," IEEE Transactions on Geoscience and Remote Sensing, submitted.
%
%If you need another two datasets (PaviaU and Salinas), please download them from http://www.ehu.eus/ccwintco/index.php/Hyperspectral_Remote_Sensing_Scenes
%
%PaviaU: http://www.ehu.eus/ccwintco/uploads/e/ee/PaviaU.mat, http://www.ehu.eus/ccwintco/uploads/5/50/PaviaU_gt.mat
%Salinas: http://www.ehu.eus/ccwintco/uploads/a/a3/Salinas_corrected.mat, http://www.ehu.eus/ccwintco/uploads/f/fa/Salinas_gt.mat


clc;clear all; 
close all;

addpath(genpath('utilities\.'));
addpath('.\data');

dataSetName='PaviaU';%PaviaU  Indianpines
switch dataSetName
case 'Indianpines'
    load Indian_pines_corrected;load Indian_pines_gt;
    data = indian_pines_corrected;
    gth  = indian_pines_gt;
case 'PaviaU'
    load PaviaU;load PaviaU_gt;
    data = paviaU;
    gth  = paviaU_gt;
end
%%=====Classifier selection====%%
classfy_options=[];
%      classfy_options.method = 'SVM';
  classfy_options.method = 'KNN';
 classfy_options.knn_k=1;%knn分类器k的大小
 classfy_options.distancekey='Distance';
 classfy_options.distancevalue='cosine';

%%==========parameter setting==========%%
demen_options=[];
demen_options.lamuda = 0.01;        % optimal parameters for lamuda =0.0001 on Indianpines and lamuda =0.01 on paviaU
demen_options.beta =0.1;            % optimal parameters for beta =0.01 on Indianpines and beta =0.1 on paviaU
demen_options.ratio =0.1;           % optimal parameters for ratio=0 on Indianpines and ratio =0.1 on paviaU
nTrEachClass = 30 ;                 %number of training data from each class
maxDim =30;                         %dimension of the feature-reduced space
radius=7;                           %Filter size
for nSeed = 1:1                     %randomly choose the training data
    [ train_x,clean_train_x,train_y,test_x,test_y]= ChooseRSdata_zxx(data,gth,nTrEachClass,nSeed,radius);
    train_x = train_x';
    clean_train_x = clean_train_x';
    test_x = test_x';
    train_y = train_y';
    test_y = test_y';
    [d_train_x,d_test_x] = RNSRP(train_x,clean_train_x,test_x,demen_options,maxDim); 
    [correctRate,optC,optG,predictY] = classifymethod(d_train_x',train_y',d_test_x',test_y',classfy_options);
    fprintf('DS:%s train:%d test:%d DR:RNRP lamuda :%f beta:%f   to %d accuracy is %f\n',dataSetName,size(train_x,2),size(test_x,2),demen_options.lamuda,demen_options.beta,maxDim,correctRate,optC,optG); 
end
