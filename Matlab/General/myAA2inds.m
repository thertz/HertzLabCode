
function [indMat] = myAA2inds(aminoAcidMat,aminoOrder)

indMat = zeros(size(aminoAcidMat));

N = length(aminoOrder);

for i=1:N
  aaInds = find(aminoAcidMat == aminoOrder{i});
  indMat(aaInds) = i;
end


if(~isempty(find(indMat == 0)))
  error('matrix does not contain amino acids, exiting..');
end