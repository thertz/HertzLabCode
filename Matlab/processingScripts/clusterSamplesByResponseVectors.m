
function [distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(respMat,numClusters,minThreshold,dist,linkageType)

% clusters responses using the response matrix across the set of selected probes using linkage. defaults to complete
% linkage with Spearman based distances. numClusters must be provided.

if(~exist('minThreshold','var'))
  minThreshold = 1000;
end

if(~exist('dist','var'))
  dist = 'spearman';
end

if(~exist('linkageType','var'))
  linkageType = 'complete';
end

%now threshold responses for correlation purposes to get rid of all noise....
maxResp       = max(respMat);
respInds      = find(maxResp > minThreshold); % this is an arbitrary threshold at this point.
distVec       = pdist(respMat(:,respInds),dist);
Zstruct       = linkage(distVec,linkageType);
distMat       = squareform(distVec);
clusterLabels = cluster(Zstruct,'maxClust',numClusters);

