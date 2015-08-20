
function [figInd] = plotIndividualResponsesForSingleAntigen(respVec,antigenName,treatmentLabels,treatmentGroups,pepData,figInd,fontSize,groupColors)

if(~exist('groupColors','var'))
  groupColors = {'blue','red'};
end

inds1 = find(treatmentLabels == 0);
inds2 = find(treatmentLabels == 1);

%handle cases where some samples are missing by inserting zeros!
n1 = length(inds1);
n2 = length(inds2);
keyboard;
if(n1 > n2)
  respMat = [respVec(inds1)' [zeros(n1-n2,1);respVec(inds2)']];
elseif(n1 < n2)
  respMat = [respVec(inds1)' zeros(n2-n1,1); respVec(inds2)'];
end


figure(figInd);
figInd = figInd + 1;
keyboard;
figHandle = bar(respMat,'grouped');

l = max(n1,n2)*2;
for i=1:l;
  if(i <= l/2)
    set(figHandle(i),'faceColor',groupColors{1});
    set(figHandle(i),'edgeColor','black');
  else
    set(figHandle(i),'faceColor',groupColors{2});
    set(figHandle(i),'edgeColor','black');
  end
end

a = gca;

%set(a,'XtickLabel',sigInds)
set(a,'FontSize',fontSize);

set(a,'YLim',[10,16]);
ylabel('log_2 MFI');

title(['Antigen specific responses for antigen',antigenName,' by groups']);

