%new function that can parse out the HA arrays printed in Scripps for the McGargill 
%antigens were spotted. Old H5N1 peptides were also spotted on this array. works on multiple arrays.
function [seqGPRstruct,otherGPRstruct] = processHIV_Scott_GPRstruct_Multiplex(gprStruct,arrayNumber,numArrays)

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
%start of zach
pepNames = seqGPRstruct.Names;
controlInds = strmatch('Dy',pepNames); % all control spots with dylight
HivInds  =  strmatch('pep',pepNames); % all the peptides/proteins

%now extract the relevanat gpr indices into the two main structs
controlGPRstruct = seqGPRstruct;
controlGPRstruct.Data    = controlGPRstruct.Data(controlInds,:); 
controlGPRstruct.Blocks  = controlGPRstruct.Blocks(controlInds); 
controlGPRstruct.Columns = controlGPRstruct.Columns(controlInds);
controlGPRstruct.Rows    = controlGPRstruct.Rows(controlInds);
controlGPRstruct.Names   = controlGPRstruct.Names(controlInds);
controlGPRstruct.IDs     = controlGPRstruct.IDs(controlInds);


seqGPRstruct.Data    = seqGPRstruct.Data(HivInds,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(HivInds,:);
seqGPRstruct.Columns = seqGPRstruct.Columns(HivInds,:);
seqGPRstruct.Rows    = seqGPRstruct.Rows(HivInds,:);
seqGPRstruct.Names   = seqGPRstruct.Names(HivInds,:);
seqGPRstruct.IDs     = seqGPRstruct.IDs(HivInds,:);




HivPepNames = seqGPRstruct.Names;
%now convert all dengue peptide names into numerical values to sort:
for i=1:length(HivPepNames)
  pepNums(i) = str2num(HivPepNames{i}(5:end));
end

% sort by peptide numbers
[Y1 I1] = sort(pepNums);
seqGPRstruct.Data    = seqGPRstruct.Data(I1,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(I1); 
seqGPRstruct.Columns = seqGPRstruct.Columns(I1);
seqGPRstruct.Rows    = seqGPRstruct.Rows(I1);
seqGPRstruct.Names   = seqGPRstruct.Names(I1);
seqGPRstruct.IDs     = seqGPRstruct.IDs(I1);










otherGPRstruct = []; % for consistency


