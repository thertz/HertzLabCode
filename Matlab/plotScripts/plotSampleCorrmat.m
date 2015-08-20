
function [figInd] = plotSampleCorrmat(figInd,corrMat,labels,colorRange,fontSize)

if(~exist('colorRange','var'))
  colorRange = [-1 1];
end

figure(figInd);
figInd = figInd+1;

imagesc(corrMat,colorRange),colorbar; 
a = gca;
set(a,'FontSize',fontSize);
set(a,'Xtick',[1:length(labels)]);
set(a,'Ytick',[1:length(labels)]);

labels = strrep(labels,'_',' ');
set(a,'XtickLabel',labels);
xticklabel_rotate([],45);
set(a,'YtickLabel',labels);
