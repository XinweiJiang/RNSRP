function [ correctRate,optC,optG, retYYSvm] = classifymethod(train_x,train_y,test_x,test_y,classfy_options )

addpath('./Tools');
if ~isfield(classfy_options,'method')
    error('Please set a classifier');
end
switch upper(classfy_options.method)
    case 'KNN'
        knn_k = classfy_options.knn_k;
        distancekey = classfy_options.distancekey;
        distancevalue = classfy_options.distancevalue;
        mdl =ClassificationKNN.fit(train_x,train_y,distancekey,distancevalue,'NumNeighbors',knn_k); 
        retYYSvm = predict(mdl,test_x);
        correctRate=length(find(retYYSvm-test_y==0))/length(test_y);
        optC=0;
        optG=0;
    case 'SVM'
        addpath('./svm2');             
        [accSvm, retYYSvm, optC, optG] = svmc(train_x, train_y, test_x, test_y, 0, 0, 1, 1);
        correctRate = accSvm;

end
end

