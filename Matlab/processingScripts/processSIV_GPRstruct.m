
function [seqGPRstruct,otherGPRstruct] = processSIV_GPRstruct(gprStruct)

nonSeqInds = [];
for i=1:length(gprStruct.Names)
  if(~isempty(findstr(gprStruct.Names{i},':')))
    nonSeqInds = [ nonSeqInds i];
  end
end
nonSeqInds = [nonSeqInds strmatch('Cy3',gprStruct.Names)'];

seqInds = setdiff([1:length(gprStruct.Names)],nonSeqInds);

seqGPRstruct = gprStruct;
seqGPRstruct.Data(nonSeqInds,:)  = [];
seqGPRstruct.Blocks(nonSeqInds)  = [];
seqGPRstruct.Columns(nonSeqInds) = [];
seqGPRstruct.Rows(nonSeqInds)    = [];
seqGPRstruct.Names(nonSeqInds)   = [];
seqGPRstruct.IDs(nonSeqInds)     = [];

otherGPRstruct = gprStruct;
otherGPRstruct.Data(seqInds,:)  = [];
otherGPRstruct.Blocks(seqInds)  = [];
otherGPRstruct.Columns(seqInds) = [];
otherGPRstruct.Rows(seqInds)    = [];
otherGPRstruct.Names(seqInds)   = [];
otherGPRstruct.IDs(seqInds)     = [];


pepNames = seqGPRstruct.Names;
otherNames = otherGPRstruct.Names;
for i=1:length(pepNames)
  pepNums(i) = str2num(pepNames{i});
end
[Y1 I1] = sort(pepNums);

seqGPRstruct.Data    = seqGPRstruct.Data(I1,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(I1); 
seqGPRstruct.Columns = seqGPRstruct.Columns(I1);
seqGPRstruct.Rows    = seqGPRstruct.Rows(I1);
seqGPRstruct.Names   = seqGPRstruct.Names(I1);
seqGPRstruct.IDs     = seqGPRstruct.IDs(I1);

[Y2 I2] = sort(otherNames);

otherGPRstruct.Data    = otherGPRstruct.Data(I2,:); 
otherGPRstruct.Blocks  = otherGPRstruct.Blocks(I2);
otherGPRstruct.Columns = otherGPRstruct.Columns(I2);
otherGPRstruct.Rows    = otherGPRstruct.Rows(I2);
otherGPRstruct.Names   = otherGPRstruct.Names(I2);
otherGPRstruct.IDs     = otherGPRstruct.IDs(I2);

