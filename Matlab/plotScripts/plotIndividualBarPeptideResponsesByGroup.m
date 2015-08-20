
% allows plotting bar plots of responses to specific antigens using a set of experimental groups, here hand-selected.
% groupNames   - names of group for each sample
% expGroups    - name of the specific groups we would like to plot (from gropuNames),

function [figInd] = plotIndividualBarPeptideResponsesByGroup(respMat,groupNames,expGroups,antigenNames,yLims,figInd,fontSize)

numGroups = length(expGroups);
for i=1:numGroups;
  inds{i} = strmatch(expGroups{i},groupNames,'exact');
  l(i) = length(inds{i});
end

maxL = max(l); % normalize all groups to have same size so bar plots are of comparable size

figure(figInd);
figInd = figInd + 1;
for i=1:numGroups
  paddingSize = maxL - length(inds{i});
  currMat = [respMat(inds{i},:); zeros(paddingSize,size(respMat,2))];

  disp(sprintf('mean responses for group %s are: %.1f %.1f %.1f %.1f',expGroups{i},mean(respMat(inds{i},:))));

  subplot(numGroups,1,i);
  figHandle= bar(currMat');

  a = gca;
  set(a,'FontSize',fontSize);

  set(a,'XtickLabel',antigenNames)
  xlabel('peptide start position');

  set(a,'yLim',yLims);  
  ylabel('MFI');

  title(sprintf('Individual responses %s (n = %d)',strrep(expGroups{i},'_',' '),length(inds{i})));
end

set(gcf,'Position',[100,200,1800,1000]);


