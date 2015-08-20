%new function that can parse out the HA arrays printed in Scripps for the McGargill 
%antigens were spotted.
function [seqGPRstruct,controlGPRstruct] = processFLU_HA_RapaGPRstruct(gprStruct)

nonSeqInds    = [strmatch('dye',gprStruct.Names)];

seqInds = setdiff([1:length(gprStruct.Names)],nonSeqInds);

seqGPRstruct = gprStruct;
seqGPRstruct.Data(nonSeqInds,:)  = [];
seqGPRstruct.Blocks(nonSeqInds)  = [];
seqGPRstruct.Columns(nonSeqInds) = [];
seqGPRstruct.Rows(nonSeqInds)    = [];
seqGPRstruct.Names(nonSeqInds)   = [];
seqGPRstruct.IDs(nonSeqInds)     = [];

controlGPRstruct = gprStruct;
controlGPRstruct.Data(seqInds,:)  = [];
controlGPRstruct.Blocks(seqInds)  = [];
controlGPRstruct.Columns(seqInds) = [];
controlGPRstruct.Rows(seqInds)    = [];
controlGPRstruct.Names(seqInds)   = [];
controlGPRstruct.IDs(seqInds)     = [];

pepNames = seqGPRstruct.Names;
for i=1:length(pepNames)
  if(pepNames{i}(1) == '#')
    pepNums(i) = str2num(pepNames{i}(2:end));
  else
    pepNums(i) = str2num(pepNames{i});
  end
end
[Y1 I1] = sort(pepNums);

seqGPRstruct.Data    = seqGPRstruct.Data(I1,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(I1); 
seqGPRstruct.Columns = seqGPRstruct.Columns(I1);
seqGPRstruct.Rows    = seqGPRstruct.Rows(I1);
seqGPRstruct.Names   = seqGPRstruct.Names(I1);
seqGPRstruct.IDs     = seqGPRstruct.IDs(I1);


