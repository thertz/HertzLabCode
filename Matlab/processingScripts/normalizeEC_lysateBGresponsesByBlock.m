
function [arrayData] = normalizeEC_lysateBGresponsesByBlock(arrayData, EC_labels, EC_antigens)

% Removes E. coli lysate responses from all samples in an arrayData struct that were made in 
% the IVTT E coli expresssion system. Uses control spots of each block to subtract local EC lystate background.
% Control spots for each block were averaged over replicates using mean/median as specified in the processing file.
% Therefore, here normalization only subtracts the mean/median bg from each mean/median antigen within that block.
% If more than one dilution, or type of bg CONTROL appears in each block, allows to subtract these separately.

% arrayData   - standard array data struct used for antigen array data
% EC_labels   - labels of the EC_lysate only control spots on the array which can be multiple names to acoomodate 
%               concentrations - assumes this is a prefix if indices are added
% EC_antigens - list of antigen names that were made in the E. Coli system. Default is all antigens 

% returns:
  % arrayData - struct with additional field named responseMatrix_EC_subtracted in which responses of all 
  % EC antigens are subtracted by the sample specific EC_lysate control responses.


if(~exist('EC_antigens','var'))
  EC_antigenInds = setdiff([1:length(arrayData.antigenNames)]);
end

blockInds = unique(arrayData.otherBlockNumbers);

for i=1:length(arrayData)  
  for j=1:length(arrayData(i).colorFlags)
    for l=1:length(EC_labels)
      respMat = arrayData(i).responseMatrix{j};    
      for k=blockInds
        blockBgInds      = find(arrayData.otherBlockNumbers == k);
        blockAntigenInds = find(arrayData.blockNumbers == k);
       
        % There should be a SINGLE index for each negative control label - 
        % since average/median across replicates was already taken during slide processing
        EC_lysateInd = strmatch(EC_labels{l},arrayData.otherAntigenNames(blockBgInds));
        EC_responses{j}{l}(:,k) = arrayData(i).otherResponseMatrix{j}(:,EC_lysateInd);            

        % subtract the background from all antigens that are within the specific block - blockAntigenInds
        respMat(:,blockAntigenInds) = respMat(:,blockAntigenInds) - repmat(arrayData(i).otherResponseMatrix{j}(:,EC_lysateInd),1,length(blockAntigenInds)); 
      end
      arrayData(i).EC_responses = EC_responses;
      arrayData(i).responseMatrix_EC_subtracted{j}{l} = respMat;  
    end
  end
end




