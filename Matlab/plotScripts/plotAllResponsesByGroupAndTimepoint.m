
function [figInd,figHandles] = plotAllResponsesByGroupAndTimepoint(mPeptResponses,groupNames,ptids,expGroups,timePoints,antigenLabel,figInd,yLims,fontSize)

if(~exist('yLims','var'))
  yLims = [-500 50000];
end

if(~exist('fontSize','var'))
  fontSize = 12;
end


numGroups = length(expGroups);

for i=1:length(timePoints)

  figHandles(i) = figure(figInd+i-1);
  
  %all samples matching current timepoint
  inds = strmatch(timePoints{i},ptids);
  
  for j=1:numGroups

    inds1 = strmatch(expGroups{j},groupNames);
    
    currInds = intersect(inds,inds1);
    if(isempty(currInds)) % for negative controls
      disp('using all samples, no timepoint for this group')
      currInds = inds1;
  end
  
    subplot(numGroups,1,j)  
    plot(mPeptResponses(currInds,:)')
  
    a = gca;
    set(a,'FontSize',fontSize);
    set(a,'ylim',yLims);
    title([strrep(expGroups{j},'x31_vn1203_',''),antigenLabel,' responses ',timePoints{i}]);
    xlabel('peptide #');
    ylabel('MFI');
  end
end

figInd = figInd + length(timePoints);