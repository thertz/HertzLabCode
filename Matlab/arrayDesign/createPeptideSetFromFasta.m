
%creates a set of peptides with a given overlap for a set of homologus
%strains, making sure that all variants of each of the strains are covered.

function [peptideStruct] = createPeptideSetFromFasta(fastaFilename,pepLength,pepOverlap,insertionPositions,rootPath,sortMethod,outputFilename)

if(~exist('rootPath','var'))
  rootPath = '/Users/thertz/';
end

if(~exist('pepLength','var'))
  pepLength = 15;
end

if(~exist('pepOverlap','var'))
  pepOverlap = 11;
end

if(~exist('sortMethod'))
  sortMethod = 'none';
end

if(~exist('outputFilename','var'))
  sepInd = strfind(fastaFilename,'.');
  outputFilename = [fastaFilename(1:sepInd-1),'_peptides_L',num2str(pepLength),'_O',num2str(pepOverlap),'.txt'];
end

[Seqs] = fastaread([rootPath,fastaFilename]);

pepGap = pepLength - pepOverlap;
currInd = 1;
for s=1:length(Seqs)
  
  seqL = length(Seqs(s).Sequence);
  begInds = [1:pepGap:seqL-pepOverlap+1];
  
  %handle insertions to correct for overlap. For each insertion start the beginning indices one step backwards to
  %allow for matches on all other positions to be retained (prevents large unrequired peptide sets due to insertions
  if(~isempty(insertionPositions{s}))
   
    for l=insertionPositions{s}
      cInd = (l+pepLength-1);
      begInds(begInds > cInd) = begInds(begInds > cInd) +1;
    end
  end
  
  for i=1:length(begInds)
  
    currBegInd = begInds(i);
    endInd     = min((begInds(i)+pepLength-1),seqL);

    peptideStruct(currInd).Header = Seqs(s).Header;
    peptideStruct(currInd).begInd   = currBegInd;
    peptideStruct(currInd).endInd   = endInd;
  
    peptideStruct(currInd).sequence = Seqs(s).Sequence(currBegInd:endInd);
    currInd = currInd + 1;
  end   
end
  
% indexStrains = [1:length(Seqs)];  
% for j=1:length(indexStrains)
%     currInds = strmatch(Seqs(indexStrains(j)).Header,{peptideStruct.Header});
    
%     if(j==1) 
%       inds = currInds;
%     else
%       currPeptides = {peptideStruct(currInds).sequence};
%       [C I] = setdiff(currPeptides,peptides);
%       inds = [inds; currInds(sort(I))]; 
%     end
%     peptides = {peptideStruct(inds).sequence};
    
% end
% peptideStruct = peptideStruct(inds);  


switch(sortMethod)

 %sort by begInd
 case 'begInd'
  [Y I] = sort([peptideStruct.begInd]);
  peptideStruct = peptideStruct(I);
  
  %sort by header and then begInd
  case 'header' 
   [Y I] = sort({peptideStruct.Header})
   peptideStruct = peptideStruct(I);

 case 'none'
end


%now write to tab delimited file:

fd = fopen([rootPath,outputFilename],'w');
fprintf(fd,'%s\t%s\t%s\t%s\t\n','Strain','beg Ind','end Ind','sequence');

for i=1:length(peptideStruct)
  fprintf(fd,'%s\t%d\t%d\t%s\t\n',peptideStruct(i).Header,peptideStruct(i).begInd,peptideStruct(i).endInd,peptideStruct(i).sequence);
end

fclose(fd);
  
  





