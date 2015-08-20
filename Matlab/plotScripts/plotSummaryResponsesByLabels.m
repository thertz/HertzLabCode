
% Plot summarized responses by group - where summary can be mean/median of each group. Based on plotAllResponsesByGroup.m
% Here - respones are grouped by labels - a numerical vector in which each group is labeled by a different integer (not restricted to specific values)

function [figInd] = plotSummaryResponsesByGroup(mPeptResponses,expGroupLabels,expGroupNames,antigenLabel,method,figInd,yLims,fontSize,lineColors)

if(~exist('yLims','var'))
  yLims = [-500 60000];
end

if(~exist('fontSize','var'))
  fontSize = 12;
end

if(~exist('method','var'))
  method = 'mean';
end

if(~exist('colors','var'))
  lineColors = ['rbkmgcy']; % assumes no more than 7 colors
end

groupLabels = unique(expGroupLabels);

figure(figInd);
figInd = figInd+1;

for i=1:length(groupLabels)

  subplot(length(groupLabels),1,i)
  inds = find(expGroupLabels == groupLabels(i))
  disp(sprintf('found %d samples from group %s',length(inds),expGroups{i}));
  
  switch(method)
    case 'mean'
      plot(mean(mPeptResponses(inds,:)),lineColors(i));
    case 'median'
      plot(median(mPeptResponses(inds,:)),lineColors(i));
  end

  a = gca;
  set(a,'FontSize',fontSize);
  set(a,'ylim',yLims);
  set(a,'yTick',[10:16])
  set(a,'yTickLabel',[10:16]);
  grid on;

  title([sprintf('%s  %s responses (n = %d)',strrep(expGroupNames{i},'_',' '),antigenLabel,length(inds))]);
  xlabel('Antigen #');
  ylabel('MFI');
  
end

set(gcf,'Position',[100,200,1800,1000]);
