
function [seqGPRstruct,otherGPRstruct] = processFLUhaGPRstruct(gprStruct)

seqInds    = strmatch('#',gprStruct.Names);
nonSeqInds = setdiff([1:length(gprStruct.Names)],seqInds);

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


pepNames = strrep(seqGPRstruct.Names,'#','');

% for i=1:length(pepNames)
%   sepInds = strfind(pepNames{i},'_');
%   pepNames{i} = pepNames{i}(1:sepInds(1)-1);
% end

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

%now further sort by dilutions:
% pepSet = unique(pepNums);
% for i=1:length(pepSet)
%   inds = strmatch(['#',num2str(i),'_'],seqGPRstruct.Names);
%   [Y I] = sort(seqGPRstruct.Names(inds));

%   seqGPRstruct.Data(inds,:)    = seqGPRstruct.Data(inds(I),:); 
%   seqGPRstruct.Blocks(inds)  = seqGPRstruct.Blocks(inds(I)); 
%   seqGPRstruct.Columns(inds) = seqGPRstruct.Columns(inds(I));
%   seqGPRstruct.Rows(inds)    = seqGPRstruct.Rows(inds(I));
%   seqGPRstruct.Names(inds)   = seqGPRstruct.Names(inds(I));
%   seqGPRstruct.IDs(inds)     = seqGPRstruct.IDs(inds(I));
  
% end
  
otherNames = strrep(otherGPRstruct.Names,'#','');
[Y2 I2] = sort(otherNames);

otherGPRstruct.Data    = otherGPRstruct.Data(I2,:); 
otherGPRstruct.Blocks  = otherGPRstruct.Blocks(I2);
otherGPRstruct.Columns = otherGPRstruct.Columns(I2);
otherGPRstruct.Rows    = otherGPRstruct.Rows(I2);
otherGPRstruct.Names   = otherGPRstruct.Names(I2);
otherGPRstruct.IDs     = otherGPRstruct.IDs(I2);

