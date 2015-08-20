
function [figInd] = plotAllResponsesByGroup(mPeptResponses,groupNames,expGroups,antigenLabel,figInd,yLims,fontSize)

if(~exist('yLims','var'))
  yLims = [-500 60000];
end

if(~exist('fontSize','var'))
  fontSize = 12;
end

numGroups = length(expGroups);

figure(figInd);
figInd = figInd+1;

for i=1:numGroups

  inds = strmatch(expGroups{i},groupNames);
  disp(sprintf('found %d samples from group %s',length(inds),expGroups{i}));
  %disp(inds);
  subplot(numGroups,1,i)  
  plot(mPeptResponses(inds,:)')
  
  a = gca;
  set(a,'FontSize',fontSize);
  set(a,'ylim',yLims);

  %set(a,'yTick',[10:16])
  %set(a,'yTickLabel',[10:16]);
  
  title([sprintf('%s  %s responses (n = %d)',strrep(expGroups{i},'_',' '),antigenLabel,length(inds))]);
  xlabel('Antigen #');
  ylabel('MFI');
  grid on;
end

set(gcf,'Position',[100,200,1800,1000]);


