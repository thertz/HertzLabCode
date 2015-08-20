
%compute the PSSM scores of the data - which are assumed to be peptides
%from the PSSM space using the given profile and symbol table

function [scores] = profileScore(data,profile,symbolTable)

if(iscell(data))
    data = cell2mat(data);
end

%convert symbolTable to amino acids the right way:
numericData = zeros(size(data));
for i=1:length(symbolTable)
  inds = find(data == symbolTable(i));
  numericData(inds) = i;
end

resMat = [];
for i=1:size(data,2)
    resMat(:,i) = profile(numericData(:,i),i);
end

scores = -sum(resMat,2);

