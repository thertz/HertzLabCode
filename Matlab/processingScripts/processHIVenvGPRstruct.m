
function [seqGPRstruct,otherGPRstruct] = processHIVenvGPRstruct(gprStruct)

seqInds    = strmatch('#',gprStruct.Names);
if(isempty(seqInds))
  seqInds = strmatch('8',gprStruct.Names);
end
nonSeqInds = setdiff([1:length(gprStruct.Names)],seqInds);

%junk labels that need to be removed from Gal file
badVec = strfind(gprStruct.Names,':');
badInds = [];
for i=1:length(badVec)
  if(~isempty(badVec{i}))
    badInds = [badInds i];
  end
end

tempInds = [seqInds ; badInds'];

seqGPRstruct = gprStruct;
seqGPRstruct.Data(nonSeqInds,:)  = [];
seqGPRstruct.Blocks(nonSeqInds)  = [];
seqGPRstruct.Columns(nonSeqInds) = [];
seqGPRstruct.Rows(nonSeqInds)    = [];
seqGPRstruct.Names(nonSeqInds)   = [];
seqGPRstruct.IDs(nonSeqInds)     = [];

otherGPRstruct = gprStruct;
otherGPRstruct.Data(tempInds,:)  = [];
otherGPRstruct.Blocks(tempInds)  = [];
otherGPRstruct.Columns(tempInds) = [];
otherGPRstruct.Rows(tempInds)    = [];
otherGPRstruct.Names(tempInds)   = [];
otherGPRstruct.IDs(tempInds)     = [];

pepNames = strrep(seqGPRstruct.Names,'#','');
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


otherNames = strrep(otherGPRstruct.Names,'#','');
[Y2 I2] = sort(otherNames);

otherGPRstruct.Data    = otherGPRstruct.Data(I2,:); 
otherGPRstruct.Blocks  = otherGPRstruct.Blocks(I2);
otherGPRstruct.Columns = otherGPRstruct.Columns(I2);
otherGPRstruct.Rows    = otherGPRstruct.Rows(I2);
otherGPRstruct.Names   = otherGPRstruct.Names(I2);
otherGPRstruct.IDs     = otherGPRstruct.IDs(I2);

