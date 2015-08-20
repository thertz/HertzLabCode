
%run the ADT prediction method
function [bindingEnergies] = ADTpredict(HLAFilename,seqFilename,fastaFlag,pepLength,outputFilename,ADTPath)

if(exist('ADTPath','var'))
  addpath(genpath(ADTPath))
end

if(~exist('pepLength','var'))
  pepLength = 9;
end


if(fastaFlag)
  Seq = fastaread(seqFilename);
else  
  fd = fopen(seqFilename,'r');
  C = textscan(fd,'%s');
  Seq.Sequence = cell2mat(C{1});
  fclose(fd);
end

if(strcmp(HLAFilename(end-3:end),'.mat'))
  load(HLAFilename);
  HLAs = alleleList;
else

  fd = fopen(HLAFilename,'r');
  [C] = textscan(fd,'%s\n');
  fclose(fd);

  adtHLAs = strrep(C{1},'*','_');
  HLAs    = C{1};
end

for i=1:length(Seq)
  [bindingEnergies{i}] = computeThreadingModelEnergies(HLAs,Seq(i).Sequence,pepLength);
end

if(exist('outputFilename','var'))
  fopen(outputFilename,'w')
  for i=1:size(bindingEnergies,2)
    fprintf(fd,'%s\t',HLAs{i});
  end
  fprintf(fd,'\n');

  for i=1:size(bindingEnergies,1)
      for j=1:size(bindingEnergies,2)
        fprintf(fd,'%.3f\t',bindingEnergies(i,j));
      end
      fprintf(fd,'\n');
  end
  fclose(fd);
end
