       
function [figInd] = plotAllResponsesByPtid(mPeptResponses,ptids,currPtid,projectName,figInd,yLims,fontSize)
%plots responses to all antigens of the specific ptid name provided. Can also plot multipled ptids
%
% mPeptResponses - responses from array as outputed by the arrayData struct
% ptids - list of ptids from arrayData struct
% currPtid - specific set of ptids that are to be plotted - cell array

if(~exist('yLims','var'))
  yLims = [-500 50000];
end

if(~exist('fontSize','var'))
  fontSize = 12;
end

figure(figInd);
figInd = figInd+1;


for i=1:length(currPtid)
  inds = strmatch(currPtid,ptids);% removed {} from after currPtid to make work
end

plot(mPeptResponses(inds,:)')

a = gca;
set(a,'FontSize',fontSize);
set(a,'ylim',yLims);
xlabel('peptide #');
ylabel('MFI');
legend(strrep(ptids(inds),'_',' '))
grid on;



