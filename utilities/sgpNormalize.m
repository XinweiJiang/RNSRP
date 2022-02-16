function [ xx ] = sgpNormalize( x)
            origBias = mean(x, 1);
            origScale = 1./sqrt(var(x, 1));
            xx = x - repmat(origBias, size(x, 1), 1);
            xx = xx.*repmat(origScale, size(xx, 1), 1);







