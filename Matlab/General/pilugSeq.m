
function [table]= pilugSeq(data)
%pilug function over amino acids (or sequence column)

if(isempty(data))
  table=[];
  return
end

minVal= min(data);
maxVal= max(data);

temp= unique(data);
if(size(temp,1)==1)  
  table(:,1)= temp';
else
  table(:,1)= temp;
end

for i=1:length(temp)
  currCount = num2str(length(find(data==table(i,1))));
  currL = length(currCount)
  begInd = 5;
  table(i,[begInd:begInd+currL-1])= currCount;
end
    