
% computes mean or median responses of fg-bg (and can easily support more) for
% a given color of interst. 
% averages by names of antigens, not assuming specific number of replicates or ordering. Simply searches by name, so
% struct need not be ordered (although it typically is).
% modified Dec. 17th 2014 to also return block number of each antigen 
% (assumes each antigen is printed in a single block at this point)
function [meanResponses,antigenNames,blockNumbers] = computeMeanResponsesByNames(gprStruct,typeFlag,colorFlag)

if(~exist('colorFlag','var'))
  colorFlag = '635';
end

if(~exist('typeFlag','var'))
  typeFlag = 'mean';
end

%extract by names: and retain order.
[antigenNames I J] = unique(gprStruct.Names);
[Y I1] = sort(I);
antigenNames = antigenNames(I1);

numPepts = length(antigenNames);

%spots are averaged, its always what to do with replicates that changes.
meanFGs = magetfield(gprStruct,['F',colorFlag,'_Mean']);
meanBGs = magetfield(gprStruct,['B',colorFlag,'_Mean']);

for i=1:numPepts
  inds = strmatch(antigenNames{i},gprStruct.Names,'exact'); 
  blockNumbers(i) = gprStruct.Blocks(inds(1)); 
  switch(typeFlag)
    case 'mean'   
     meanResponses(i) = mean( meanFGs(inds) - meanBGs(inds) );
   case 'median'
    meanResponses(i) = median( meanFGs(inds) - meanBGs(inds) );
  end
end
