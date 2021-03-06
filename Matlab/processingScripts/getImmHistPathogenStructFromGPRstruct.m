
function [pathogenStruct] = getImmHistPathogenStructFromGPRstruct(seqGPRstruct,numRepeats)

%computes infor for of each pathogen on a seqGPRstruct for Exp1, assuming pathogens have been
%ordered. Provides info on where the pathogen probes are, their names etc.
% assumes treplicates were used for all samples/
%
% numRepeats - number of replicates for each reagent on slide (typically 3)

pathogenNames = unique(seqGPRstruct.IDs); 


for i=1:length(pathogenNames)
  
  currInds = strmatch(pathogenNames{i},seqGPRstruct.IDs,'exact');
  
  pathogenStruct(i).name = pathogenNames{i};
  pathogenStruct(i).inds = currInds;
  pathogenStruct(i).begInd = ceil(min(currInds)/numRepeats); %beginning index of pathogen after averaging over treplicates.
  pathogenStruct(i).endInd = max(currInds)./numRepeats;
  
  for j=1:length(currInds)
    pathogenStruct(i).antigens(j).name = seqGPRstruct.Names{currInds(j)};
    repInds = strmatch(pathogenStruct(i).antigens(j).name,seqGPRstruct.Names,'exact');
    
    %sanity check - verify that each antigen is really repaeted numRepeat times!
    %if(length(repInds) ~= numRepeats)
    %  error('found an antigen with less than expected number of repeats in GPR struct, exiting');
    %end
    pathogenStruct(i).antigens(i).ind  = repInds;
  end
end

%now resort by begin indices on seqGPRstruct
[Y I] = sort([pathogenStruct.begInd]);
pathogenStruct = pathogenStruct(I);