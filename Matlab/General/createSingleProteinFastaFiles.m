
%creates a sepearte fasta file for each of the proteins from the given fasta file.

function [] = createSingleProteinFastaFiles(fastaFile,saveDir,headerSep,sepInd)

if(~exist('headerSep','var'))
  headerSep = '|';
  sepInd    = 1;
end

[data] = fastaread(fastaFile);

for i=1:length(data)
  sepInds = findstr(data(i).Header,headerSep);
  
  if(sepInd == 1)
    protName = data(i).Header(1:sepInds(1)-1);
  else
    protName = data(i).Header(sepInds(sepInd-1)+1:sepInds(sepInd)-1);
  end
    
  fd = fopen([saveDir,protName,'.fasta'],'w');
  fprintf(fd,'>%s\n%s\n',data(i).Header,data(i).Sequence);
  fclose(fd);
end
  