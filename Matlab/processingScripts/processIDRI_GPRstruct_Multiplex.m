%designed for array/gpr printed on 3-11-14
function [seqGPRstruct,otherGPRstruct] = processIDRI_GPRstruct_Multiplex(gprStruct,arrayNumber,numArrays)

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


%sort by strain and then by pep number (According to indexed list)
strainNames = seqGPRstruct.Names;
pathogenNames = {'Pep'};% names of strain in array, add other names only if placed on the same gal file,

pathogenInds = [];  
for i=1:length(pathogenNames)
  
  inds = strmatch(pathogenNames{i},strainNames);
  strainNums = [];
  for j=1:length(inds)
    strainNums(j) = str2num(strainNames{inds(j)}(length(pathogenNames{i})+2:end)); %numbers start in second spot after the above pathogen names 
  end
  [Y1 I1] = sort(strainNums);
  
  pathogenInds = [pathogenInds; inds(I1)];
end


seqGPRstruct.Data    = seqGPRstruct.Data(pathogenInds,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(pathogenInds); 
seqGPRstruct.Columns = seqGPRstruct.Columns(pathogenInds);
seqGPRstruct.Rows    = seqGPRstruct.Rows(pathogenInds);
seqGPRstruct.Names   = seqGPRstruct.Names(pathogenInds);
seqGPRstruct.IDs     = seqGPRstruct.IDs(pathogenInds);

otherGPRstruct = []; % for consistency




