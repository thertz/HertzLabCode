function [Energies] = get_sequence_HLA_binding_energies_by_PSSMs(Seq,PSSMs,peptide_len)

if(size(Seq,1) == 1)
  [peptides]=get_peptide_matrix_from_sequence(Seq,peptide_len);
else
  peptides = Seq;
end
[Energies]=get_gral_HLA_peptides_energies_by_PSSM(peptides,PSSMs);

