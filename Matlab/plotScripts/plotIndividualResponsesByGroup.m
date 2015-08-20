
function [figInd,sigInds] = plotIndividualResponsesByGroup(respMat,treatmentLabels,treatmentGroups,antigen_pValues,pThreshold,strainInds,...
  strainName,pepData,figInd,fontSize,groupColors)

if(~exist('groupColors','var'))
  groupColors = {'blue','red'};
end

sigInds    = intersect(find(antigen_pValues < pThreshold),strainInds);

respMat = respMat(:,sigInds);

inds1 = find(treatmentLabels == 0);
inds2 = find(treatmentLabels == 1);
%handle cases where some samples are missing by inserting zeros!
n1 = length(inds1);
n2 = length(inds2);

if(n1 > n2)
  respMat = [respMat ;zeros(n1-n2,size(respMat,2))]';
elseif(n1 < n2)
  respMat = [respMat(inds1,:) ;zeros(n2-n1,size(respMat,2)); respMat(inds2,:)]';
end

figure(figInd);
figInd = figInd + 1;
figHandle= bar(respMat);

l= size(respMat,2);
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
set(a,'Xtick',[1:length(sigInds)]);
set(a,'XtickLabel',[pepData(sigInds).begInd])

%set(a,'XtickLabel',sigInds)
set(a,'FontSize',fontSize);
xlabel('peptide start position');

%set(a,'YLim',[10,16]);
%ylabel('log_2 MFI');
ylabel('MFI');

title(['Antigen specific responses by groups to strain ',strainName]);

