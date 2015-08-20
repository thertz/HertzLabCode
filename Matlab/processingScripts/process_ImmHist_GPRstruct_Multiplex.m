%process_ImmHist_GPRstruct_Multiplex  is made from processFlu......
%new function that can parse out the HA arrays printed in Scripps for the McGargill 
%antigens were spotted. Old H5N1 peptides were also spotted on this array. works on multiple arrays.
function [seqGPRstruct,otherGPRstruct] = process_ImmHist_GPRstruct_Multiplex(gprStruct,arrayNumber,numArrays)

numBlocks = max(gprStruct.Blocks);
numBlocksPerArray = numBlocks/numArrays;

begBlock = (arrayNumber-1)*numBlocksPerArray+1;
endBlock = begBlock + numBlocksPerArray-1;

begInds = find(gprStruct.Blocks == begBlock);
endInds = find(gprStruct.Blocks == endBlock);

seqGPRstruct = gprStruct;
seqGPRstruct.Data    = seqGPRstruct.Data([begInds(1):endInds(end)],:);
seqGPRstruct.Blocks  = seqGPRstruct.Blocks([begInds(1):endInds(end)]);
seqGPRstruct.Columns = seqGPRstruct.Columns([begInds(1):endInds(end)]);
seqGPRstruct.Rows    = seqGPRstruct.Rows([begInds(1):endInds(end)]);
seqGPRstruct.Names   = seqGPRstruct.Names([begInds(1):endInds(end)]);
seqGPRstruct.IDs     = seqGPRstruct.IDs([begInds(1):endInds(end)]);

% the "dye" spots will be given a peptide number of -1, so can be easily removed below
pepNames = seqGPRstruct.Names;
for i=1:length(pepNames)
  if(strcmp(pepNames{i},'dye'))
    pepNums(i) = -1;
  else
    pepNums(i) = str2num(pepNames{i});
  end
end
[Y1 I1] = sort(pepNums);
%now remove all dye inds
dyeInds = find(Y1==-1);
I1(dyeInds) = [];

seqGPRstruct.Data    = seqGPRstruct.Data(I1,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(I1); 
seqGPRstruct.Columns = seqGPRstruct.Columns(I1);
seqGPRstruct.Rows    = seqGPRstruct.Rows(I1);
seqGPRstruct.Names   = seqGPRstruct.Names(I1);
seqGPRstruct.IDs     = seqGPRstruct.IDs(I1);

otherGPRstruct = []; % for consistency


