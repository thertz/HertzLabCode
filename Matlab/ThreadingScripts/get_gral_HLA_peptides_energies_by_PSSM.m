function [Energies]=get_gral_HLA_peptides_energies_by_PSSM(peptides,PSSMs);

num_aa=20;
[num_peptides,peptide_len] = size(peptides);

if ischar(peptides)
   peptides = get_aa_indexes(peptides);    
end

num_HLAs = length(PSSMs);
Energies=zeros(num_peptides,num_HLAs);

offset = (0:peptide_len-1)*(num_aa+1);

peptides = peptides+ones(num_peptides,1)*offset;

for i=1:num_HLAs
    Energies(:,i)=sum(PSSMs{i}(peptides),2);
end

Energies=-(Energies+PSSMs{1}(end,end));