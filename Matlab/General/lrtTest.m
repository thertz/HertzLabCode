
function [LRT] = lrtTest(countMat)

LRT = 0;
for i=1:2
  for j=1:2
    
    n = countMat(i,j)./sum(countMat(i,:));
    m = sum(countMat(:,j))./sum(countMat(:));
    LRT = LRT + countMat(i,j)*log(n/m);
  end
end
LRT = 2*LRT;
