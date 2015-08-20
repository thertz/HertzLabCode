
function [antigenFilters] = computeRespMatLabelBlindedFilters(respMat,filterNames,filterThresholds,posThresholds)

%create label blinded filters that allow only analyzing a subset of antigens on the array for group comparisons to
%allow handling multiplicity

for i=1:length(filterNames)
  
  switch(filterNames{i})
    
   case 'PercentResponders' 
    
    fMat = (respMat > posThresholds(i));
    pVec = sum(fMat,1)./size(fMat,1);
    antigenFilters(i,:) = pVec > filterThresholds{i};
    
  end
end