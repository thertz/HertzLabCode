
function [figInd] = plotAllResponsesByClusters(mPeptResponses,clusterLabels,antigenLabel,figInd,yLims,fontSize)

if(~exist('yLims','var'))
  yLims = [-500 60000];
end

if(~exist('fontSize','var'))
  fontSize = 12;
end

uniqLabels = unique(clusterLabels);
numGroups = length(uniqLabels);

figure(figInd);
figInd = figInd+1;

xLims = [0 size(mPeptResponses,2)];

for i=1:numGroups

  inds = find(clusterLabels == uniqLabels(i));
  disp(sprintf('found %d samples from cluster %d',length(inds),i));
 
  subplot(numGroups,1,i)  
  plot(mPeptResponses(inds,:)')
  
  a = gca;
  set(a,'FontSize',fontSize);
  set(a,'ylim',yLims);
  set(a,'xLim',xLims)

  %set(a,'yTick',[10:16])
  %set(a,'yTickLabel',[10:16]);
  
  title([sprintf('cluster %d  %s responses (n = %d)',i,antigenLabel,length(inds))]);
  %xlabel('Antigen #');
  ylabel('MFI');
  grid on;
end

set(gcf,'Position',[100,200,1800,1000]);


