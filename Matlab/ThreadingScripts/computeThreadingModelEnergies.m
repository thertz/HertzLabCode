% runs the Threading model for a given set of alleles and provides binding
% energies for every peptide of length pepLength along the given protein
% sequence.
% 
% input:
% alleleList - a list of alleles. alleles can be double digit or four digits. notation using '_' instead of '*'
% example: alleleList = {'A_0201' 'B_2705'};  
% 
% protein - the given protein sequence
% pepLength (optional) - the length of the peptides considered.
% added a flag to allow different prediction engines for A's and B alleles.
% Modified 7.2.08, protein can now be a matrix of peptides of a fixed
% size. However in this case does not hanlde non AA chars like '-' OR 'X'
function [bindingEnergies] = computeThreadingModelEnergies(alleleList,protein,pepLength,useLocusSpecificEngineFlag,...
						  mutateAlleleFlag,mutPosition,mutAminoAcid,HLAsStruct,HLAsNames,bindingEngine)

if(~exist('pepLength','var'))
    pepLength = 9;
end

if(~exist('mutateAlleleFlag','var'))
  mutateAlleleFlag = 0;
end

if(~exist('useLocusSpecificEngineFlag','var'))
  useLocusSpecificEngineFlag = 0;
end

% load the Threading model trained by Manuel
%load DataFiles/exp_binding_energies_benchmark9_assymetric_no_out  
%NEWER MODEL - not using LANL data???
if(~exist('bindingEngine','var'))
  if(useLocusSpecificEngineFlag)
    if(~isempty(strmatch('A',alleleList)))
      load DataFiles\exp_binding_eng_Bench_HLA_A;
      bindingEngine = exp_binding_eng_Bench_HLA_A;
    else
      load DataFiles\exp_binding_eng_Bench_HLA_B;
      bindingEngine = exp_binding_eng_Bench_HLA_B;
    end
  else
    load DataFiles/exp_binding_eng_Bench_no_out_all_cluster_1
    bindingEngine = exp_binding_eng_Bench_no_out_all_cluster_1;
    
    %load DataFiles/exp_binding_eng_Bench_no_out_all_cluster_1_three_or_less;
    %bindingEngine = exp_binding_eng_Bench_no_out_all_cluster_1_three_or_less;
    
    %load DataFiles/exp_binding_energies_benchmark9_assymetric_no_out;
    %bindingEngine = exp_binding_energies_benchmark9_assymetric_no_out;

  end
end


%if this is a protein and not a list of peptides
badPepInds = [];
if(size(protein,1) == 1)
  
  %identify empty letters in sequence.
  badInds = strfind(protein,'-');

  for i=1:length(badInds)
    badPepInds = [badPepInds badInds-8:badInds];
  end
  
  %now temporarily stick some letter into the bad positions!
  protein = strrep(protein,'-','A');

  % transform protein into all pepLength peptides.
  [protPeptides]=get_peptide_matrix_from_sequence(protein,pepLength); 
 
else
  protPeptides = protein;
end

if(~exist('HLAsStruct','var')) 
  %loads  All_HLA_names  and HLA_Strings. Loading outside the function
  %load \\manuelrg1\ProteinFoldingProject\DataFiles\HLA_Strings_struct;
  load HLA_Strings_struct;
  HLAsStruct = HLA_Strings;
  HLAsNames  = All_HLA_names;
end

alleleInds = find_HLA_indexes(HLAsNames,alleleList);

if(mutateAlleleFlag)

  %compute "Energies" for all HLAs on "HLA_Strings(indexes)" changing aminoacids on positions in "pos" by aminoacids on "aas"
  [bindingEnergies] = get_gral_HLA_peptides_energies_by_struct_sust_aas(protPeptides,HLAsStruct(alleleInds),...
										     bindingEngine,mutPosition,mutAminoAcid);
else
  [bindingEnergies] = get_gral_HLA_peptides_energies_by_struct(protPeptides,HLAsStruct(alleleInds),bindingEngine);
end

bindingEnergies(badPepInds,:) = ones(length(badPepInds),1)*mean(bindingEnergies); %some default value for the time being!
