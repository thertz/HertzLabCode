
function [figHandle] = plotAntigenSpecificPvalues(groupResponses,type,threshold,figInd)

if(~exist('threshold','var'))
  threshold = 0.05;
end

if(~exist('figInd','var'))
  figInd = 1;
end


numPlots = length(groupResponses)

figure(figInd);
figInd = figInd+1;

for i=1:length(groupResponses)
  
  subplot(numPlots,1,i);
  hold on;
  
  switch(type)
   case 'ranksum'
    pValues        = groupResponses(i).antigenSpecific_pValuesRankSum;
   case 'Fisher'
    pValues = groupResponses(i).antigenSpecific_pValuesFisher;
  end
   
  antigenFilters = groupResponses(i).antigenFilters;  
  inds = find(antigenFilters);
  plot(inds,-log10(pValues(inds)),'o');
  
  line([inds; inds],[zeros(1,length(inds));-log10(pValues(inds))]);
  
  plot([1:length(antigenFilters)],-log10(threshold),'k-');
  
  a = gca;
  set(a,'FontSize',14);
  set(a,'ytick',[-log10(1) -log(0.5) -log10(0.1) -log10(0.05) -log10(0.01) -log10(0.001) -log10(0.0001)]);
  set(a,'ytickLabel', [1 0.5 0.1 0.05 0.01 0.001 0.0001]);
  %set(a,'Ylim',[1 -log10(0.0001)]);
      
  title(groupResponses(i).name);
  
end

