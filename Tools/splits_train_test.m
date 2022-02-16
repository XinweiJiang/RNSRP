function [train_x,train_y,test_x,test_y] = splits_train_test(dataset,oldLabels,options)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
%     dataset   each column is a data point;
%     oldlabels  row ,each clumn is a label;
%     trainratio select train size 
%     selecttype random normal 
%     oldLabels=[1;1;1;1;1;1;2;2;1;3;2;1;2;2;3;3;3];
%     trainratio = 0.5;
    if ~exist('options', 'var')
        options.selecttype = 'random'; 
    end
    unique_labelN = unique(oldLabels);
    label_index = [1:size(oldLabels, 2)];
    train_index = [];
    test_index = [];
    switch options.selecttype
        case 'random'
            for i = 1:length(unique_labelN)  %遍历所有的标签
                %找到当前标签为对应的的序号
                index_i = find(oldLabels == unique_labelN(i));
                if isempty(index_i) % 如果该类别在数据集中不存在
                    continue;
                end
                %随机函数
                rp = randperm(length(index_i)); % random permutation
                %floor是取整函数
                rp_ratio = rp(1:floor(length(index_i)*options.trainratio));
                %提取当前对应的序号的数据出来
                ind = index_i(rp_ratio);
                train_index = [train_index,ind];
        %         Ctrain = [Ctrain,length(ind)];
        %         Ctest = [Ctest,length(setdiff(index_i,ind))];
                mm=unique_labelN(i);
                fprintf('label %d train_num:%d .total:%d\n',mm,size(ind,2),size(index_i,2));
            end
        case 'normal'
            for i = 1:length(unique_labelN)  %遍历所有的标签
                %找到当前标签为对应的的序号
                index_i = find(oldLabels == unique_labelN(i));
                if isempty(index_i) % 如果该类别在数据集中不存在
                    continue;
                end
                %floor是取整函数
                %提取当前对应的序号的数据出来
                ind = index_i(options.selectnum);
                train_index = [train_index,ind];
        %         Ctrain = [Ctrain,length(ind)];
        %         Ctest = [Ctest,length(setdiff(index_i,ind))];
                mm=unique_labelN(i);
                fprintf('label %d train_num:%d .total:%d\n',mm,size(ind,2),size(index_i,2));
            end
        otherwise
            error('Wrong paramaters');
    end
    test_index = [test_index;setdiff(label_index,train_index)];
    fprintf('totaltrain_num:%d totaltest_num:%d \n',size(train_index,2),size(test_index,2));
    
    train_x = dataset(:,train_index);
    train_y = oldLabels(:,train_index); 
    test_x = dataset(:,test_index);
    test_y = oldLabels(:,test_index);
    
%     mean_x = mean(train_x,2);
%     train_num = size(train_x,2);
%     test_num = size(test_x,2);
%     train_x = train_x - repmat(mean_x,1,train_num);
%     test_x = test_x - repmat(mean_x,1,test_num);
end

