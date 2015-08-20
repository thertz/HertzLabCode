%new function that can parse out the Dengue arrays printed in Scripps for Eva Harris
%antigens were spotted. Currently used for single arrays, but can be used for more.
% seqGPRstruct - all dengue peptides from all 4 serotypes
% controlGPRstruct - all control peptides that are Ab contact sites and their muntants.

function [seqGPRstruct,controlGPRstruct] = processDENGUE_GPRstruct_Multiplex(gprStruct,arrayNumber,numArrays)

numBlocks = max(gprStruct.Blocks);
numBlocksPerArray = numBlocks/numArrays;

begBlock = (arrayNumber-1)*numBlocksPerArray+1;
endBlock = begBlock + numBlocksPerArray-1;

begInds = find(gprStruct.Blocks == begBlock);
endInds = find(gprStruct.Blocks == endBlock);
keyboard;
seqGPRstruct = gprStruct;
seqGPRstruct.Data    = seqGPRstruct.Data([begInds(1):endInds(end)],:);
seqGPRstruct.Blocks  = seqGPRstruct.Blocks([begInds(1):endInds(end)]);
seqGPRstruct.Columns = seqGPRstruct.Columns([begInds(1):endInds(end)]);
seqGPRstruct.Rows    = seqGPRstruct.Rows([begInds(1):endInds(end)]);
seqGPRstruct.Names   = seqGPRstruct.Names([begInds(1):endInds(end)]);
seqGPRstruct.IDs     = seqGPRstruct.IDs([begInds(1):endInds(end)]);


% the "HA" and "FLAG" spots will be given a peptide number of -1, so can be easily removed below
pepNames = seqGPRstruct.Names;
controlInds = strmatch('Control',pepNames); % all control peptides start with Control
dengueInds  = strmatch('D',pepNames); % all Dengue peptides start with D

%now extract the relevanat gpr indices into the two main structs
controlGPRstruct = seqGPRstruct;
controlGPRstruct.Data    = controlGPRstruct.Data(controlInds,:); 
controlGPRstruct.Blocks  = controlGPRstruct.Blocks(controlInds); 
controlGPRstruct.Columns = controlGPRstruct.Columns(controlInds);
controlGPRstruct.Rows    = controlGPRstruct.Rows(controlInds);
controlGPRstruct.Names   = controlGPRstruct.Names(controlInds);
controlGPRstruct.IDs     = controlGPRstruct.IDs(controlInds);


seqGPRstruct.Data    = seqGPRstruct.Data(dengueInds,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(dengueInds); 
seqGPRstruct.Columns = seqGPRstruct.Columns(dengueInds);
seqGPRstruct.Rows    = seqGPRstruct.Rows(dengueInds);
seqGPRstruct.Names   = seqGPRstruct.Names(dengueInds);
seqGPRstruct.IDs     = seqGPRstruct.IDs(dengueInds);


denguePepNames = seqGPRstruct.Names;
%now convert all dengue peptide names into numerical values to sort:
for i=1:length(denguePepNames)
  pepNums(i) = str2num(denguePepNames{i}(2:end));
end
% sort by peptide numbers
[Y1 I1] = sort(pepNums);
seqGPRstruct.Data    = seqGPRstruct.Data(I1,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(I1); 
seqGPRstruct.Columns = seqGPRstruct.Columns(I1);
seqGPRstruct.Rows    = seqGPRstruct.Rows(I1);
seqGPRstruct.Names   = seqGPRstruct.Names(I1);
seqGPRstruct.IDs     = seqGPRstruct.IDs(I1);


%order the controlGPRstruct:
cInds = [strmatch('Control1',controlGPRstruct.Names) ;strmatch('Control2',controlGPRstruct.Names)];
controlGPRstruct.Data    = controlGPRstruct.Data(cInds,:); 
controlGPRstruct.Blocks  = controlGPRstruct.Blocks(cInds); 
controlGPRstruct.Columns = controlGPRstruct.Columns(cInds);
controlGPRstruct.Rows    = controlGPRstruct.Rows(cInds);
controlGPRstruct.Names   = controlGPRstruct.Names(cInds);
controlGPRstruct.IDs     = controlGPRstruct.IDs(cInds);



