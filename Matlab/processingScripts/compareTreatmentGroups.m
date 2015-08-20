
function [groupStats] = compareTreatmentGroups(respMat,labels,posThreshold)

%defines which responses are considered below threshold
if(~exist('posThreshold','var'))
	posThreshold = log2(1024);
end

% Compare overall responses:

[pValue,totResponses] = compareGroupsByTotalArrayResponses(respMat,labels);
groupStats.totResponses_pValue  = pValue;
groupStats.totResponses         = totResponses;


%now compare responses of each antigen separately:

[pValuesRankSum,pValuesFisher,group1responseHigher,fisherCounts] = compareGroupsByAntigen(respMat,labels,posThreshold);
		
groupStats.singleAntigen_pValuesRankSum        = pValuesRankSum;
groupStats.singleAntigen_pValuesFisher         = pValuesFisher;
groupStats.FisherCounts                        = fisherCounts; 
groupStats.singleAntigens_group1responseHigher = group1responseHigher;

