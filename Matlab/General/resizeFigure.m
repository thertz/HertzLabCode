
function [figHandle] = resizeFigure(figHandle,xSize,ySize)

if(~exist('xSize','var'))
  xSize =  24;
end
  
if(~exist('ySize','var'))
  ySize = 18;
end

set(gcf,'PaperUnits','centimeters');
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;

set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[100 100 xSize*50 ySize*50]);

   