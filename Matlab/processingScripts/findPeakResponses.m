
function [peakResponses] = findPeakResponses(responseMatrix,threshold)

% returns a binay array with all responses that are above a certain threshold for each row in the responseMatrix data

if(~exist('threshold','var'))
  threshold = 20000;
end

peakResponses = responseMatrix > threshold;
