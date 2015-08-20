% process_ImmHist_GPRstruct_Multiplex  is made from processFl
% new function that can parse multiple arrays
% note that pathogen list is hard coded in this file at this point -
% currently includes HSV_1,HSV_2,VZV and KSHV
function [seqGPRstruct,otherGPRstruct] = processImmHistVer1_GPRstruct_Multiplex(gprStruct,arrayNumber,numArrays)

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


%sort by pathogen and then by orf number (According to indexed list)
protNames = seqGPRstruct.Names;
pathogenNames = {'HSV_1','HSV_2','VZV','KSHV','Neg'};

pathogenInds = [];
for i=1:length(pathogenNames)
  
  inds = strmatch(pathogenNames{i},protNames);
  protNums = [];
  for j=1:length(inds)
    protNums(j) = str2num(protNames{inds(j)}(length(pathogenNames{i})+2:end));
  end
  [Y1 I1] = sort(protNums);
  
  pathogenInds = [pathogenInds; inds(I1)];
end

seqGPRstruct.Data    = seqGPRstruct.Data(pathogenInds,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(pathogenInds); 
seqGPRstruct.Columns = seqGPRstruct.Columns(pathogenInds);
seqGPRstruct.Rows    = seqGPRstruct.Rows(pathogenInds);
seqGPRstruct.Names   = seqGPRstruct.Names(pathogenInds);
seqGPRstruct.IDs     = seqGPRstruct.IDs(pathogenInds);

otherGPRstruct = []; % for consistency


