
function [figHandle] = printFig(figHandle,figName,xSize,ySize,jpgFlag)

% Prints figures to a file in two formats: eps and png (both lossless)
% Params:
% figHandle - handle to figure (figure's index)
% figName   - full name of figure including full path
% xSize, ySize (optional) - the size (width, height) of the figure

figure(figHandle);

if(~exist('xSize','var'))
  xSize =  24;
end
  
if(~exist('ySize','var'))
  ySize = 18;
end

if(~exist('jpgFlag','var'))
	pdfFlag = 0;
end

set(gcf,'PaperUnits','centimeters');
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;

set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[100 100 xSize*50 ySize*50]);

if(jpgFlag)
	print('-djpeg90',[figName,'.jpg']);
else
	print('-depsc2','-painters',[figName,'.eps']);
	print('-dpng','-r600',[figName,'.png']);
end
