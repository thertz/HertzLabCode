
function [figInd] = plotResponseDendrogram(Zstruct,labels,figInd,fontSize,bwFlag,titleStr)

if(~exist('fontSize','var'))
  fontSize = 12;
end

if(~exist('bwFlag','var'))
  bwFlag = 0;
end

if(~exist('titleStr','var'))
  titleStr = [];
end


figure(figInd)
figInd = figInd + 1;
h = dendrogram(Zstruct,0,'Orientation','right','labels',labels);

if(bwFlag)
  for i=1:length(h)
    set(h,'Color','black');
  end
end
a = gca;
set(a,'FontSize',fontSize);

if(~isempty(titleStr))
  title(titleStr);
end

set(gcf,'Position',[100,200,1400,1000]);