
%removes the renowned asterix char from all fasta files within the given
%directory

function [] = removeAsterixFromFastaFiles(fastaDir)

files = dir([fastaDir,'*.fasta']);

for i=1:length(files)
  [data] = fastaread([fastaDir,files(i).name]);
  
  astInds = strfind(data.Sequence,'*');
  if(length(astInds) > 1)
    disp(sprintf('more than one asterix removed from file %s',files(i).name));
  end
  data.Sequence(astInds) =[];
  
  fd = fopen([fastaDir,files(i).name],'w');
  fprintf(fd,'>%s\n',data.Header);
  fprintf(fd,'%s\n',data.Sequence);
  fclose(fd);
end
