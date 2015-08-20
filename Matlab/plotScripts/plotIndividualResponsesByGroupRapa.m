
function [figInd,figHandle,sigInds] = plotIndividualResponsesByGroupRapa(groupResponses,pValue,strainInds,strainName,origArrayPepData,figInd,fontSize,groupColors)

if(~exist('groupColors','var'))
  groupColors = {'blue','red'};
end

sigInds    = intersect(find(groupResponses.antigenSpecific_pValuesFisher < pValue),strainInds);

respMat1 = groupResponses.responseMat{1}(:,sigInds);
respMat2 = groupResponses.responseMat{2}(:,sigInds);

%handle cases where some samples are missing by inserting zeros!
n1 = size(respMat1,1);
n2 = size(respMat2,1);

if(n1==n2)
  respMat = [respMat1;respMat2]';
elseif(n1>n2)
  respMat = [respMat1; [respMat2;zeros(n1-n2,size(respMat2,2))]]';
else
  respMat = [[respMat1;zeros(n2-n1,size(respMat1,2))]; respMat2]';
end

figure(figInd);
figInd = figInd + 1;
figHandle= bar(respMat);

for i=1:n1+n2
  if(i <= n1)
    set(figHandle(i),'faceColor',groupColors{1});
    set(figHandle(i),'edgeColor','black');
  else
    set(figHandle(i),'faceColor',groupColors{2});
    set(figHandle(i),'edgeColor','black');
  end
end

a = gca;
set(a,'Xtick',[1:length(sigInds)]);
set(a,'XtickLabel',[origArrayPepData(sigInds).begInd])

%set(a,'XtickLabel',sigInds)
set(a,'FontSize',fontSize);
xlabel('peptide start position');
ylabel('MFI');
%title(['Antigen specific responses by groups to strain ',strainName]);

groupNames = [groupResponses(1).groupNames(1:2)];
groupNames = strrep(groupNames,'x31_vn1203_','');

%legend(figHandle([1,n1+1]),groupNames);
