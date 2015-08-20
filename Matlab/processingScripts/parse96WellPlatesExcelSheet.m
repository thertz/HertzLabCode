

%read excell sheet with 96 well plates in standard format and parse them out
%to allow creating 384 well plates. Assumes fixed format where title of plate
%is always: Plate X - ...... and plate begins exactly 2 lines after (title row).

function [plateMaps] = parse96WellPlatesExcelSheet(plateFilename,rootPath)

 
if(~exist('rootPath','var'))
  rootPath =  '/Users/thertz/';
end

%[NUMERIC,TXT,RAW]=xlsread([rootPath,plateFilename]);

fd = fopen([rootPath,plateFilename],'r');
i=1;
while(~feof(fd))
  C{i} = fgetl(fd);
  i = i+1;
end
fclose(fd);

%find all titles of plates - how many plates are there
inds = strmatch('Plate ',C);

for i=1:length(inds)
  plateMaps(i) = parsePlateMap(C(inds(i):inds(i)+10));
end
