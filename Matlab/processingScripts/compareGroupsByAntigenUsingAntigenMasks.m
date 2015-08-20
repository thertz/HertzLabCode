
function [antigenStats] = compareGroupsByAntigenUsingAntigenMasks(respMat,dataLabels,antigenMasksStruct,filterNames,filterThresholds,minResponseThreshold)

% Filter antigens using a pre-defined filter on blinded data - currently % responders
% Compute Fisher's exact p-values on responders vs. non-responders for each antigen separately
% Compute q-Values for p-values, limiting analysis to antigens that belong to the filtered antigen list
% Returns antigenStats- a struct with various stat measures from this procedure.

for i=1:length(antigenMasksStruct)

  % compute label blinded filtering of antigens across the array - currenlty using %responders - set at 80%
  cmdStr = ['[antigenFilter] = computeRespMatLabelBlindedFilters(respMat(:,antigenMasksStruct(i).inds),filterNames,filterThresholds,minResponseThreshold);'];
  eval(cmdStr)
  
  antigenStats(i).name = antigenMasksStruct(i).name;
  antigenStats(i).filteredAntigenInds = find(antigenFilter); %list of indices to antigens in which more than % repsonderes responded ( blinded)

  cmdStr = ['[stats] = compareTreatmentGroups(respMat(:,antigenMasksStruct(i).inds),dataLabels,minResponseThreshold);'];
  
  eval(cmdStr);
  antigenStats(i).totResponses_pValue                 = stats.totResponses_pValue;
  antigenStats(i).totResponses                        = stats.totResponses;
  antigenStats(i).singleAntigen_pValuesRankSum        = stats.singleAntigen_pValuesRankSum;            
  antigenStats(i).singleAntigen_pValuesFisher         = stats.singleAntigen_pValuesFisher;
  antigenStats(i).FisherCounts                        = stats.FisherCounts;
  antigenStats(i).singleAntigens_group1responseHigher = stats.singleAntigens_group1responseHigher; 

  antigenStats(i).sigInds = intersect(antigenStats(i).filteredAntigenInds,find(antigenStats(i).singleAntigen_pValuesFisher < 0.05));
  antigenStats(i).pValues = antigenStats(i).singleAntigen_pValuesFisher(antigenStats(i).sigInds);

  antigenStats(i).sigIndsRankSum = intersect(antigenStats(i).filteredAntigenInds,find(antigenStats(i).singleAntigen_pValuesRankSum < 0.05));
  antigenStats(i).pValuesRankSum = antigenStats(i).singleAntigen_pValuesRankSum(antigenStats(i).sigIndsRankSum);

  % remove NANs from ranksum p-values created due to lack of responses
  antigenStats(i).singleAntigen_pValuesRankSum(find(isnan(antigenStats(i).singleAntigen_pValuesRankSum))) = 1;

  % Q-value correction - adjust all p-values that are within the set of filtered inds
  [FDR, qValues] = mafdr(antigenStats(i).singleAntigen_pValuesFisher(antigenStats(i).filteredAntigenInds));
  antigenStats(i).qInds   = antigenStats(i).filteredAntigenInds(qValues < 0.2);
  antigenStats(i).qValues = qValues(qValues < 0.2);

  [FDR, qValues] = mafdr(antigenStats(i).singleAntigen_pValuesRankSum(antigenStats(i).filteredAntigenInds));
  antigenStats(i).qIndsRankSum   = antigenStats(i).filteredAntigenInds(qValues < 0.2);
  antigenStats(i).qValuesRankSum = qValues(qValues < 0.2);

end
