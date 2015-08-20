
function [pValue,totResponses] = compareGroupsByTotalArrayResponses(responseMat,labels)

% compute two-sided wilcoxon ranksum test to compare the total response of two groups marked by labels 
% Responses are summed across all array probes, and are assumed to be BG subtracted.

%first take respMat and set all responses below zero to zero:
responseMat(responseMat<0)=0;

totResponses = sum(responseMat,2);

uniqLabels = unique(labels);

inds1 = find(labels == uniqLabels(1));
inds2 = find(labels == uniqLabels(2));

[pValue] = ranksum(totResponses(inds1),totResponses(inds2));
