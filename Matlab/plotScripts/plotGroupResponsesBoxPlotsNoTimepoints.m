% version using groupResponse struct that does not contain timepoints (flat format)
function [figInd,figHandle] = plotGroupResponsesBoxPlotsNoTimepoints(groupResponses,groupInds,antigenInds,AbType,figInd,fontSize)

if(~exist('fontSize','var'))
  fontSize = 12;
end

if(~exist('figInd','var'))
  figInd = 1;
end

labels  = [];
respVec = [];
for j=groupInds
  respVec = [respVec sum(groupResponses(j).responseMat(:,antigenInds),2)'];
  labels = [labels ones(1,size(groupResponses(j).responseMat,1))*j ];
end
  
figHandle = figure(figInd);
h= notBoxPlot(respVec,labels,[],'sdline');
a = gca;
set(a,'Xtick',[1:length(groupInds)]);
set(a,'XtickLabel',{groupResponses.name});
  

figInd = figInd + 1;

[pValue] = ranksum(respVec(labels==groupInds(1)),respVec(labels==groupInds(2)));
  
a = gca;
set(a,'FontSize',fontSize);
title([AbType,' p=',num2str(pValue)]);
  
ylabel('MFI');

