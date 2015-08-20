
% Plot summarized responses by unupservised clusters - where summary can be mean/median of each group. Based on plotAllResponsesByGroup.m

function [figInd] = plotSummaryResponsesByClusters(mPeptResponses,clusterLabels,antigenLabel,method,figInd,yLims,fontSize,lineColors)

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

numGroups = length(unique(clusterLabels));

figure(figInd);
figInd = figInd+1;

for i=1:numGroups

  subplot(numGroups,1,i)
  inds = find(clusterLabels == i);
  disp(sprintf('found %d samples from cluster %d',length(inds),i));

  %single response is summarized by default...  
  if(length(inds)==1)
    plot(mPeptResponses(inds,:),lineColors(i));
  else

    switch(method)
      case 'mean'
        plot(mean(mPeptResponses(inds,:)),lineColors(i));
      case 'median'
        plot(median(mPeptResponses(inds,:)),lineColors(i));
    end
  end
  a = gca;
  set(a,'FontSize',fontSize);
  set(a,'ylim',yLims);
  set(a,'yTick',[10:16])
  set(a,'yTickLabel',[10:16]);
  grid on;

  title([sprintf('cluster %d %s responses (n = %d)',i,antigenLabel,length(inds))]);
  xlabel('Antigen #');
  ylabel('MFI');
  
end

set(gcf,'Position',[100,200,1800,1000]);
