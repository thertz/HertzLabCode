
function [bgResponses] = computeBgResposnesFromArrayDataStruct(arrayData,bgLabel)

% compute spot specific responses for negative controls used to estimated background for positivity calls. Computes
% these separately for each arrayData struct entry (each entry is a different exeprimental day).
% Negative controls refers to samples here, and not control antigens.
%
% bgLabel - allows various types of negative controls to be used and analyzed separately

for i=1:length(arrayData)
  
  bgResponses(i).bgLabel = bgLabel;
  
  inds = strmatch(bgLabel,arrayData(i).ptids);
  if(~isempty(inds))
    
    for j=1:length(arrayData(i).colorFlags)
      bgResponses(i).colorFlag{j}    = arrayData(i).colorFlags{j};
      bgResponses(i).respMat{j}      = arrayData(i).responseMatrix{j}(inds,:);
      bgResponses(i).avgResp(j,:)    = mean(arrayData(i).responseMatrix{j}(inds,:));
      bgResponses(i).medianResp(j,:) = median(arrayData(i).responseMatrix{j}(inds,:));
      bgResponses(i).stdResp(j,:)    = std(arrayData(i).responseMatrix{j}(inds,:));
    end
  end
end



