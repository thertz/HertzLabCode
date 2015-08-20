
function [figInd] = barPlotGroups(mPeptResponses,groupNames,expGroups,antigenLabel,figInd,yLims, fontSize)

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
  bar(mPeptResponses(inds,:)')
  
  a = gca;
  set(a,'FontSize',fontSize);
  set(a,'ylim',yLims);
  title([strrep(expGroups{i},'_',' '),' ',antigenLabel,' responses']);
  xlabel('peptide #');
  ylabel('MFI');
  grid on;
end

set(gcf,'Position',[100,200,1800,1000]);


