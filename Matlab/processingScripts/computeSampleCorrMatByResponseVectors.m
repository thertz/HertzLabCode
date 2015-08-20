
function [corrMat] = computeSampleCorrMatByResponseVectors(respMat,minThreshold,corrType)

% computes correlation of response vectors, limiting only to probes in which the maximal response is above the
% minThreshold. Defaults to Spearman Correlation.

if(~exist('minThreshold','var'))
  minThreshold = 1000;
end

if(~exist('corrType','var'))
  corrType = 'Spearman';
end

%now threshold responses for correlation purposes to get rid of all noise....
maxResp       = max(respMat);
respInds      = find(maxResp > minThreshold); % this is an arbitrary threshold at this point.

corrMat = corr(respMat(:,respInds)','type',corrType);


