
function [arrayData] = normalizeEC_lysateBGresponses(arrayData, EC_label, method, EC_antigens)

% Removes E. coli lysate responses from all samples in an arrayData struct that were made in 
% the IVTT E coli expresssion system. 

% arrayData   - standard array data struct used for antigen array data
% EC_label    - label of the EC_lysate only control spots on the array
% method      - method used for normalization - can be 'median', 'mean' and 'max'
% EC_antigens - list of antigen names that were made in the E. Coli system. Default is all antigens 

% returns:
  % arrayData - struct with additional field named responseMatrix_EC_subtracted in which responses of all 
  % EC antigens are subtracted by the sample specific EC_lysate control responses.


EC_lysateInds = strmatch(EC_label,arrayData.antigenNames);

if(~exist('EC_antigens','var'))
  EC_antigenInds = setdiff([1:length(arrayData.antigenNames)],EC_lysateInds);
else
  EC_antigenInds = [];
  for i=1:length(EC_antigens)
    EC_antigenInds = [EC_antigenInds ; strmatch(EC_antigens{i},arrayData(1).antigenNames)];
  end
end
  
for i=1:length(arrayData)  
  for j=1:length(arrayData(i).colorFlags)

      switch(method)
        case 'median'
            EC_responses{j} = median(arrayData(i).responseMatrix{j}(:,EC_lysateInds),2);  
        case 'mean'
            EC_responses{j} = mean(arrayData(i).responseMatrix{j}(:,EC_lysateInds),2);  
        case 'max'
            EC_responses{j} = max(arrayData(i).responseMatrix{j}(:,EC_lysateInds),[],2);  
      end
      
      respMat = arrayData(i).responseMatrix{j};
      respMat(:,EC_antigenInds) = respMat(:,EC_antigenInds) - repmat(EC_responses{j}',size(EC_antigenInds))';
      arrayData(i).responseMatrix_EC_subtracted{j} = respMat;      
  end
  arrayData(i).EC_responses = EC_responses;
end




