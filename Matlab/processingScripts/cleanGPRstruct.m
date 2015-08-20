%modified to now include all the possible forms of empty cells that can be cleaned away.
function [gprStruct] = cleanGPRstruct(gprStruct)

inds = unique([strmatch('"empty"',gprStruct.IDs)' strmatch('"empty"',gprStruct.Names)' strmatch('na',gprStruct.Names)'...
	      strmatch('"unknown"',gprStruct.Names)' strmatch('"0"',gprStruct.Names)' ]);

gprStruct.Data(inds,:)  = [];
gprStruct.Blocks(inds)  = [];
gprStruct.Columns(inds) = [];
gprStruct.Rows(inds)    = [];
gprStruct.Names(inds)   = [];
gprStruct.IDs(inds)     = [];

inds = strmatch('""',gprStruct.IDs);

gprStruct.Data(inds,:)  = [];
gprStruct.Blocks(inds)  = [];
gprStruct.Columns(inds) = [];
gprStruct.Rows(inds)    = [];
gprStruct.Names(inds)   = [];
gprStruct.IDs(inds)     = [];

gprStruct.Names = strrep(gprStruct.Names,'"','');
gprStruct.IDs = strrep(gprStruct.IDs,'"','');









