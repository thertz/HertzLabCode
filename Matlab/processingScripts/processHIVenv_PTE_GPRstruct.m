%new function that deals with both conB and PTE peptides, and also with
%additional protein reagents.
function [seqGPRstruct,pteGPRstruct,otherGPRstruct] = processHIVenv_PTE_GPRstruct(gprStruct)

conBseqInds = strmatch('8',gprStruct.Names)';
%remove all PTE labels - all conB labels are 4 digit long! (in future we need
%to distinguish them with # or another sign!!!!
remInds = [];
for i=1:length(conBseqInds)
  if(length(gprStruct.Names{conBseqInds(i)}) ~= 4)
    remInds = [remInds i];
  end
end
conBseqInds(remInds) = [];
 
PTEinds = [];
for i=1:482 %the number of PTE peptides for ENV, names by their number!
  PTEinds =  [PTEinds  strmatch(num2str(i),gprStruct.Names,'exact')'];
end

nonSeqInds = setdiff([1:length(gprStruct.Names)],[conBseqInds   PTEinds]);
 
seqGPRstruct = gprStruct;
seqGPRstruct.Data    = seqGPRstruct.Data(conBseqInds,:);
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(conBseqInds); 
seqGPRstruct.Columns = seqGPRstruct.Columns(conBseqInds);
seqGPRstruct.Rows    = seqGPRstruct.Rows(conBseqInds);
seqGPRstruct.Names   = seqGPRstruct.Names(conBseqInds);
seqGPRstruct.IDs     = seqGPRstruct.IDs(conBseqInds);
 

pteGPRstruct = gprStruct;
pteGPRstruct.Data    = pteGPRstruct.Data(PTEinds,:);
pteGPRstruct.Blocks  = pteGPRstruct.Blocks(PTEinds); 
pteGPRstruct.Columns = pteGPRstruct.Columns(PTEinds);
pteGPRstruct.Rows    = pteGPRstruct.Rows(PTEinds);
pteGPRstruct.Names   = pteGPRstruct.Names(PTEinds);
pteGPRstruct.IDs     = pteGPRstruct.IDs(PTEinds);
 

otherGPRstruct = gprStruct;
otherGPRstruct.Data    = otherGPRstruct.Data(nonSeqInds,:);
otherGPRstruct.Blocks  = otherGPRstruct.Blocks(nonSeqInds); 
otherGPRstruct.Columns = otherGPRstruct.Columns(nonSeqInds);
otherGPRstruct.Rows    = otherGPRstruct.Rows(nonSeqInds);
otherGPRstruct.Names   = otherGPRstruct.Names(nonSeqInds);
otherGPRstruct.IDs     = otherGPRstruct.IDs(nonSeqInds);
 
 
%sort peptide structs by peptide numbers
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

pteNames = strrep(pteGPRstruct.Names,'#','');
for i=1:length(pteNames)
  pteNums(i) = str2num(pteNames{i});
end
[Y2 I2] = sort(pteNums);

pteGPRstruct.Data    = pteGPRstruct.Data(I2,:); 
pteGPRstruct.Blocks  = pteGPRstruct.Blocks(I2); 
pteGPRstruct.Columns = pteGPRstruct.Columns(I2);
pteGPRstruct.Rows    = pteGPRstruct.Rows(I2);
pteGPRstruct.Names   = pteGPRstruct.Names(I2);
pteGPRstruct.IDs     = pteGPRstruct.IDs(I2);


[Y3 I3] = sort(otherGPRstruct.Names);
otherGPRstruct.Data    = otherGPRstruct.Data(I3,:); 
otherGPRstruct.Blocks  = otherGPRstruct.Blocks(I3);
otherGPRstruct.Columns = otherGPRstruct.Columns(I3);
otherGPRstruct.Rows    = otherGPRstruct.Rows(I3);
otherGPRstruct.Names   = otherGPRstruct.Names(I3);
otherGPRstruct.IDs     = otherGPRstruct.IDs(I3);


