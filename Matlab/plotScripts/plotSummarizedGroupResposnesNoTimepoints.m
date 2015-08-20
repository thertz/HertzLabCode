%version tha tworks on groups where there are no timpoints within each group (flat format of groupResponses)

function [figInd] = plotSummarizedGroupResposnesNoTimepoints(groupResponses,groupInds,summaryStat,antigenInds,antigenLabel,figInd,yLims,fontSize)

numPlots = length(groupInds);
  
figure(figInd);
hold on;
for j=1:numPlots
    
  switch(summaryStat)
    
   case 'Mean'
    respVec(j,:) = mean(groupResponses(j).responseMat(:,antigenInds));
    
   case 'Median'
    respVec(j,:) = median(groupResponses(j).responseMat(:,antigenInds));
    end
  
    subplot(numPlots,1,j)
    plot(respVec(j,:));
    a = gca;
    set(a,'FontSize',fontSize);
    set(a,'ylim',yLims);
    
    groupLabel = strrep(groupResponses(j).name,'_','');
        
    title([groupLabel,antigenLabel]);
    xlabel('peptide #');
    ylabel('MFI');
end
figInd = figInd+1;
    
    
    
  
  
  