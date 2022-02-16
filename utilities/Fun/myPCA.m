function [eigvector_PCA, eigvalue_PCA] = myPCA(data, PCARatio)
%PCA	Principal Component Analysis
%
%	Usage:
%       [eigvector, eigvalue] = PCA(data, ReducedDim)
%       [eigvector, eigvalue] = PCA(data)
% 
%             Input:
%               data       - Data matrix. Each column vector is a data point.
%
%           PCARatio   - if PCARatio<=1...if PCARatio>1...
%
%             Output:
%               eigvector_PCA - Each column is an embedding function, for a new
%                           data point (row vector) x,  y = x*eigvector
%                           will be the embedding result of x.
%               eigvalue_PCA  - The sorted eigvalue of PCA eigen-problem. 
%
%   Written by Lishan Qiao
%                                                   
[nFea,nSmp] = size(data);
data=data-data*ones(nSmp)/nSmp;
if nSmp >= nFea
    ddata = data*data';
    ddata = max(ddata,ddata');
 
    [eigvector_PCA, eigvalue_PCA] = eig(ddata);
    eigvalue_PCA = diag(eigvalue_PCA);
    clear ddata;
 
    maxEigValue = max(abs(eigvalue_PCA));
    eigIdx = find(eigvalue_PCA/maxEigValue < 1e-12);
    eigvalue_PCA(eigIdx) = [];%remove the columns which match small eigvalue
    eigvector_PCA(:,eigIdx) = [];
 
    [junk, index] = sort(-eigvalue_PCA);
    eigvalue_PCA = eigvalue_PCA(index);
    eigvector_PCA = eigvector_PCA(:, index);
        
    %=======================================
    if PCARatio > 1
        idx = PCARatio;
        if idx < length(eigvalue_PCA)
            eigvalue_PCA = eigvalue_PCA(1:idx);
            eigvector_PCA = eigvector_PCA(:,1:idx);
        end
    elseif PCARatio < 1
        sumEig = sum(eigvalue_PCA);
        sumEig = sumEig*PCARatio;
        sumNow = 0;
        for idx = 1:length(eigvalue_PCA)
            sumNow = sumNow + eigvalue_PCA(idx);
            if sumNow >= sumEig
                break;
            end
        end
        eigvalue_PCA = eigvalue_PCA(1:idx);
        eigvector_PCA = eigvector_PCA(:,1:idx);
    end
    %=======================================

else %% nSmp < nFea
    ddata = data'*data;
    ddata = max(ddata,ddata');
 
    [eigvector, eigvalue_PCA] = eig(ddata);
    eigvalue_PCA = diag(eigvalue_PCA);
    clear ddata;
 
    maxEigValue = max(eigvalue_PCA);
    eigIdx = find(eigvalue_PCA/maxEigValue < 1e-12);
    eigvalue_PCA(eigIdx) = [];
    eigvector(:,eigIdx) = [];
 
    [junk, index] = sort(-eigvalue_PCA);
    eigvalue_PCA = eigvalue_PCA(index);
    eigvector = eigvector(:, index);
        
    %=======================================
    if PCARatio > 1
        idx = PCARatio;
        if idx < length(eigvalue_PCA)
            eigvalue_PCA = eigvalue_PCA(1:idx);
            eigvector = eigvector(:,1:idx);
        end
    elseif PCARatio < 1
        sumEig = sum(eigvalue_PCA);
        sumEig = sumEig*PCARatio;
        sumNow = 0;
        for idx = 1:length(eigvalue_PCA)
            sumNow = sumNow + eigvalue_PCA(idx);
            if sumNow >= sumEig
                break;
            end
        end
        eigvalue_PCA = eigvalue_PCA(1:idx);
        eigvector = eigvector(:,1:idx);
    end
    %=======================================

    eigvalue_PCA = eigvalue_PCA.^-.5;
    eigvector_PCA = (data*eigvector).*repmat(eigvalue_PCA',nFea,1);
 
    clear eigvector;
end