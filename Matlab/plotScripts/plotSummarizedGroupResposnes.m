
% Specific plotting function used in the McGargill dataset, depracated and replaced by more general summary printing file: plotSummaryResponsesByGroup.m

function [figInd] = plotSummarizedGroupResposnes(groupResponses,groupInds,summaryStat,antigenInds,antigenLabel,figInd,yLims,fontSize)

numPlots = length(groupInds);

for i=1:length(groupResponses)
  
  figure(figInd);
  hold on;
  for j=1:numPlots
    
    switch(summaryStat)
      
     case 'Mean'
      respVec(j,:) = mean(groupResponses(i).responseMat{groupInds(j)}(:,antigenInds));
  
     case 'Median'
      respVec(j,:) = median(groupResponses(i).responseMat{groupInds(j)}(:,antigenInds));
    end
  
    subplot(numPlots,1,j)
    plot(respVec(j,:));
    a = gca;
    set(a,'FontSize',fontSize);
    set(a,'ylim',yLims);
    
    groupLabel = strrep(groupResponses(i).groupNames{groupInds(j)},'x31_vn1203_','');
    timeLabel = strrep(groupResponses(i).time,'_',' ');
    
    title([groupLabel,antigenLabel,' day ',timeLabel]);
    xlabel('peptide #');
    ylabel('MFI');
  end
  figInd = figInd+1;
end
    
    
    
  
  
  