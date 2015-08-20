
function [seqGPRstruct,otherGPRstruct] = processImmHist_GPRstruct(gprStruct)

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

%currently this is a junk GPR struct that has only cy3 and undefined spots.
otherGPRstruct = gprStruct;
otherGPRstruct.Data(seqInds,:)  = [];
otherGPRstruct.Blocks(seqInds)  = [];
otherGPRstruct.Columns(seqInds) = [];
otherGPRstruct.Rows(seqInds)    = [];
otherGPRstruct.Names(seqInds)   = [];
otherGPRstruct.IDs(seqInds)     = [];


%now we want to sort by the specific pathogens and peptides we have...
pathogenPrefix = {'HSV-1-','HSV-2-','KSHV-','EBV','CMV','FLU','H1','H3','Adeno5-','ConB-','8','E.coli_lysate','E.coli lysate'};
pathogenNames  = {'HSV-1','HSV-2','KSHV','EBV','CMV','Influenza-A','Influenza-A','Influenza-A','Adeno-5','HIV', ...
		  'HIV','E.coli_lysate','E.coli lysate'};

%sort struct by pathogen Names:
inds = [];
IDs  = {};
for i=1:length(pathogenPrefix)
  currInds  = strmatch(pathogenPrefix{i},seqGPRstruct.Names);
  currNames = seqGPRstruct.Names(currInds);
  %handle case where some pathogens were not printed
  if(isempty(currInds))
    continue;
  end
  %for HIV conB peptides, we want to sort by peptide numbers
  if(strmatch('8',pathogenPrefix{i})) %HIV CON-b names!
    for j=1:length(currNames)
      pepNums(j) = str2num(currNames{j});
    end
    [Y I] = sort(pepNums);
  else
    [Y I] = sort(currNames);
  end
  inds = [inds currInds(I)'];
  IDs  = [ IDs; repmat(pathogenNames(i),length(currNames),1) ];
end


seqGPRstruct.Data    = seqGPRstruct.Data(inds,:); 
seqGPRstruct.Blocks  = seqGPRstruct.Blocks(inds); 
seqGPRstruct.Columns = seqGPRstruct.Columns(inds);
seqGPRstruct.Rows    = seqGPRstruct.Rows(inds);
seqGPRstruct.Names   = seqGPRstruct.Names(inds);
seqGPRstruct.IDs     = IDs; %replace IDs that are tyipcally empty with the real Pathogen IDs...


