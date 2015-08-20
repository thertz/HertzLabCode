
function [figInd] = plotGroupResponsesBoxPlots(respVec,groupLabels,groupNames,pValue,figInd,fontSize,titleStr)


if(~exist('fontSize','var'))
  fontSize = 12;
end

if(~exist('figInd','var'))
  figInd = 1;
end

if(~exist('titleStr','var'))
	titleStr = '';
end

figure(figInd);
figInd = figInd + 1;

h= notBoxPlot(respVec,groupLabels,[],'sdline');
a = gca;
set(a,'Xtick',[0:length(unique(groupLabels))-1]);
set(a,'XtickLabel',groupNames);
   
a = gca;
set(a,'FontSize',fontSize);
title([titleStr,' p=',num2str(pValue)]);
  
%ylabel('log_2 MFI');
ylabel('MFI');

