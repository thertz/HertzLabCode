% process_ImmHist_GPRstruct_Multiplex  is made from processFl
% new function that can parse multiple arrays
% localBgGPRstruct contains null reactions in each block separately ordered by block
function [seqGPRstruct,localBgGPRstruct] = processKSHVOrf_GPRstruct_Multiplex(gprStruct,arrayNumber,numArrays)

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


%sort by pathogen and then by dilution)
protNames = seqGPRstruct.Names;

pathogenNames = {'K', 'DK'}; 
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

bgNames       = {'null','Dnull'}; %names of null antigens

bgInds = [];
for i=1:length(bgNames)
  
  inds = strmatch(bgNames{i},protNames);
  protNums = [];
  for j=1:length(inds)
    protNums(j) = str2num(protNames{inds(j)}(length(bgNames{i})+2:end));
  end
  [Y1 I1] = sort(protNums);
  
  bgInds = [bgInds; inds(I1)];
end

localBgGPRstruct          = gprStruct;
localBgGPRstruct.Data     = gprStruct.Data(bgInds,:); 
localBgGPRstruct.Blocks   = gprStruct.Blocks(bgInds); 
localBgGPRstruct.Columns  = gprStruct.Columns(bgInds);
localBgGPRstruct.Rows     = gprStruct.Rows(bgInds);
localBgGPRstruct.Names    = gprStruct.Names(bgInds);
localBgGPRstruct.IDs      = gprStruct.IDs(bgInds);



