
%parse a CSV entry (or any other neatly delimited string) and return all tokens in a cell array of strings
%(originally used to parse the A and B HIV CD8 epitope list). Also eats
%spaces. Removes quotation marks if such exist around string.
function [tokens] = parseCSVentry(str,delimiter)

if(~exist('delimiter','var'))
  delimiter = ',';
end

%remove  "" marks that might encompass entire string.
str = strrep(str,'"','');

tokens = {};
%parse data from entry
while(~isempty(str))
  [T,str] = strtok(str,delimiter);
  tokens = [tokens strrep(T,' ','')];
end

