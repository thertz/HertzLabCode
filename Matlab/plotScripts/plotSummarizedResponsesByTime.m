
function [figInd] = plotSummarizedResponsesByTime(groupResponses,summaryStat,antigenInds,antigenLabel,timeInds,figInd,yLims,fontSize)

numPlots = length(groupResponses);

for i=timeInds

  figure(figInd);
  hold on;
  
  for j=1:length(groupResponses)
      
    switch(summaryStat)
      
     case 'Mean'
      respVec(j,:) = mean(groupResponses(j).responseMat{i}(:,antigenInds));
  
     case 'Median'
      respVec(j,:) = median(groupResponses(j).responseMat{i}(:,antigenInds));
    end
  
    subplot(numPlots,1,j)
    plot(respVec(j,:));
    a = gca;
    set(a,'FontSize',fontSize);
    set(a,'ylim',yLims);
    
    title([strrep(groupResponses(j).groupNames{i},'x31_vn1203_',''),antigenLabel,' day ',strrep(groupResponses(j).time,'_',' ')]);
    xlabel(['peptide #');
    ylabel('MFI') 
  end
  figInd = figInd+1;
end
    
    
    
  
  
  