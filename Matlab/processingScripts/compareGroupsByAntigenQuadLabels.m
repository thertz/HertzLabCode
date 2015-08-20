
function [pValuesRankSum,pValuesFisher,group1responseHigher,fisherCounts] = compareGroupsByAntigenQuadLabels(respMat,labels,posThreshold)

% Fisher's exact test on quad labels using a positivity threshold
% anlaysis to specific antigens in which are selected in a treatment blinded manner.
%
% posThreshold - threshold to define a positive response for binarizing outpouts

%first take respMat and set all responses below zero to zero:
respMat(respMat<0)=0;

numAntigens = size(respMat,2);

uniqLabels = unique(labels);

inds1 = find(labels == uniqLabels(1));
inds2 = find(labels == uniqLabels(2));

for i=1:numAntigens
  pValuesRankSum(i) = ranksum(respMat(inds1,i),respMat(inds2,i));
  %[bla,pValues(i)] = ttest2(respMat(inds1,i),respMat(inds2,i));

  %binarize resposnes into yes/no and compute Fisher's exact for each antigen:
  a = sum(respMat(inds1,i) < posThreshold);
  b = sum(respMat(inds1,i) >= posThreshold);
  
  c = sum(respMat(inds2,i) < posThreshold);
  d = sum(respMat(inds2,i) >= posThreshold);
  
  fisherCounts(i,:) = [a b c d];

 [Ppos,Pneg,Pboth]=Fisherextest(a,b,c,d);
 pValuesFisher(i) = Pboth;
 if(b>d) % more responders in group 1
   group1responseHigher(i) = 1;
 else
   group1responseHigher(i) = 0;
 end
end
