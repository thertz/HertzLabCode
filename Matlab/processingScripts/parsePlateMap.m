
%parses cell array with a 96 well plate map with fixed format. All entries
%are tab delimited...
% title line: Plate X....
% blank line
% column numbers - 1.... 12 
% Row number and map - A p1 p2 p3 p4....p12

function [plateMap] = parsePlateMap(plateData)

%check data
if(isempty(strmatch('Plate',plateData{1})))
  error('Plate data in wrong format, exiting');
end

colNums = textscan(plateData{3},'%*d %d %d %d %d %d %d %d %d %d %d %d %d','Delimiter','\t'); 

plateMap.title = strtrim(plateData{1});
plateMap.data = cell(8,12);
[plateMap.data{:}] = deal('');
for i=1:8
  C = textscan(plateData{i+3},'%*s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter','\t'); 
  strVec = [C{:}];
  plateMap.data(i,1:length(strVec)) = strVec;
end

%make sure that all empty cells contain an empty string....

