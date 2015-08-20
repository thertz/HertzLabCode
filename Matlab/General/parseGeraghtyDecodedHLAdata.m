
%reads an HLA data file in Geraghty format and returns it in a parsable struct.
function [hlaStruct] = parseGeraghtyDecodedHLAdata(hlaFilename)

fd = fopen(hlaFilename,'r');

%SampleID	A1	A2	B1	B2	C1	C2


[C] = textscan(fd,'%d %s %s %s %s %s %s','Delimiter','\t','HeaderLines',1);
fclose(fd);

for i=1:length(C{1})
  
  hlaStruct(i).ptid = C{1}(i);
  hlaStruct(i).rawHLAs = { C{2}{i} C{3}{i} C{4}{i} C{5}{i} C{6}{i} C{7}{i} };
  
  if(isempty(hlaStruct(i).rawHLAs{5}))
    hlaStruct(i).hlas    = { C{2}{i}(1:6) C{3}{i}(1:6) C{4}{i}(1:6) C{5}{i}(1:6) };
  else
    hlaStruct(i).hlas    = { C{2}{i}(1:6) C{3}{i}(1:6) C{4}{i}(1:6) C{5}{i}(1:6) C{6}{i}(1:6) C{7}{i}(1:6) };
  end
  
end