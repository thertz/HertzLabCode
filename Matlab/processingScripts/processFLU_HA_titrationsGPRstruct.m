%new function that can parse out the HA arrays in whcih titrations of the
%antigens were spotted.
function [seqGPRstruct,controlGPRstruct] = processFLU_HA_titrationsGPRstruct(gprStruct)

nonSeqInds    = [strmatch('Cy3 IgA',gprStruct.Names) ;strmatch('PBS',gprStruct.Names) ]; 

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


%fix labels of 1mg/ml orig to avoid confusion with the 1mg/ml label
seqGPRstruct.Names = strrep(seqGPRstruct.Names,' 1mg/ml orig',' 1 mg/ml orig ');

%there are 4 titrations: 5mg/ml, 1mg/ml 1mg/ml orig and 0.5mg/ml
titrationStrs = {' 5mg/ml',' 2.5mg/ml',' 1mg/ml',' 0.5mg/ml',' 1 mg/ml orig'};

orderVec = [];
for i=1:length(titrationStrs)
  allInds = strfind(seqGPRstruct.Names,titrationStrs{i});
  tInds = [];
  for j=1:length(allInds)
    if(~isempty(allInds{j}))
      tInds = [tInds j];
    end
  end
  
  %sort all antigens from the same titration:
  currAgs = seqGPRstruct.Names(tInds);
  currAgNames = strrep(currAgs,titrationStrs{i},''); 
  
  pepNums =  [];
  for j=1:length(currAgNames)
    pepNums(j) = str2num(currAgNames{j});
  end
  [Y1 I1] = sort(pepNums);
  
  orderVec = [orderVec tInds(I1)];
  
end

seqGPRstruct.Data    = seqGPRstruct.Data(orderVec,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(orderVec); 
seqGPRstruct.Columns = seqGPRstruct.Columns(orderVec);
seqGPRstruct.Rows    = seqGPRstruct.Rows(orderVec);
seqGPRstruct.Names   = seqGPRstruct.Names(orderVec);
seqGPRstruct.IDs     = seqGPRstruct.IDs(orderVec);
  
controlNames = strrep(controlGPRstruct.Names,'#','');
[Y2 I2] = sort(controlNames);

controlGPRstruct.Data    = controlGPRstruct.Data(I2,:); 
controlGPRstruct.Blocks  = controlGPRstruct.Blocks(I2);
controlGPRstruct.Columns = controlGPRstruct.Columns(I2);
controlGPRstruct.Rows    = controlGPRstruct.Rows(I2);
controlGPRstruct.Names   = controlGPRstruct.Names(I2);
controlGPRstruct.IDs     = controlGPRstruct.IDs(I2);

