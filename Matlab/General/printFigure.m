
function [figHandle] = printFigure(figHandle,figName,typeFlags)

% Prints figures to a file in two formats: eps and png (both lossless)
% Params:
% figHandle - handle to figure (figure's index)
% figName   - full name of figure including full path
% typeFlags (optional) - which formats to print figure to. default is jpg.

figure(figHandle);

if(~exist('typeFlags','var'))
	typeFlags = {'jpg'};
end

for i=1:length(typeFlags)
	
	if(strmatch('jpg',typeFlags))
		print('-djpeg90',[figName,'.jpg']);	
	end
	
	if(strmatch('eps',typeFlags))
		print('-depsc2','-painters',[figName,'.eps']);	
	end

	if(strmatch('png',typeFlags))
		print('-dpng','-r600',[figName,'.png']);
	end
end
