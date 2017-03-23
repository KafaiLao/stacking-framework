function ndata = normalizeSignal(data, dim)
% Normalize the signal to zero mean and unit variance
% [Input]
% data: Array with any dimensions
% dim: dimension that are going to be normalized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataMean = mean(data, dim);
ShapeDataMean = ones(1,length(size(data)));
ShapeDataMean(dim) = size(data, dim);
ndata = data - repmat(dataMean, ShapeDataMean);
dataStd = std(ndata, 0, dim);
ndata = ndata./repmat(dataStd, ShapeDataMean);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Last modified 23/03/2017 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Ka Fai Lao, University of Macau %%%%%%%%%%%%%%%%%%%%%
end