%version works on groupResponses without timepoints (flat file)
function [figInd,figHandle,sigInds] = plotIndividualResponsesByGroupNoTimepoints(groupResponses,groupInds,pValue,strainInds,strainName,origArrayPepData,figInd,fontSize)

sigInds    = intersect(find(groupResponses(1).antigenSpecific_pValuesFisher < pValue),strainInds);

respMat1 = groupResponses(groupInds(1)).responseMat(:,sigInds);
respMat2 = groupResponses(groupInds(2)).responseMat(:,sigInds);

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
    %set(figHandle(i),'faceColor',[0.75,0.75,0.75]);
    set(figHandle(i),'faceColor','red');
    set(figHandle(i),'edgeColor','black');
  else
    %set(figHandle(i),'faceColor','white');
    set(figHandle(i),'faceColor','blue');
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
title(['Antigen specific responses by groups to strain ',strainName]);

groupNames = {groupResponses(groupInds).name};
groupNames = strrep(groupNames,'_',' ');

legend(figHandle([1,n1+1]),groupNames);
