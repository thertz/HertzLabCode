
%creates a matrix of (N*M)x2 where N and M are the lengths of the index
%vectors xInd and yInd. Uses matlab tricks to efficiently create all pairs
%without iterating at all. 
% xInd and yInd are both column vectors.
% if you would like to tranform a vector of this indice responses into a
% matrix where each row of this pairMat is an entry you should code
% [resMat] = reshape(resVec,N,M); 

function [pairMat] = createAllIndicePairs(xInds,yInds);

N = length(xInds);
M = length(yInds);

xInds = repmat(xInds,1,M);
yInds = repmat(yInds,N,1);

pairMat = [xInds(:),yInds(:)];


