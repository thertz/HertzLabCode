
%parse slide to sample mapping files for automated processing of data from a given experiment. File is in 3 columns:
% slide no. sample group

function [slideNumbers,sampleNames,labels,arrayNumbers,comments] = parseMultiplexSlideMappingFile(mappingFilename,commentFlag)

if(~exist('commentFlag','var'))
  commentFlag = 0;
end



fd = fopen(mappingFilename,'r');
if(commentFlag)
  [C] = textscan(fd,'%s\t%s\t%s\t%d\t%s','HeaderLines',1,'Delimiter','\t');
else
  [C] = textscan(fd,'%s\t%s\t%s\t%d\t','HeaderLines',1);
end
fclose(fd);

slideNumbers = C{2};
sampleNames  = C{1};
labels       = C{3};
arrayNumbers = C{4};
if(commentFlag)
  comments     = C{5};
end