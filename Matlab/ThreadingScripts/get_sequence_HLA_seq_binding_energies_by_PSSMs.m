function [Energies] = get_sequence_HLA_seq_binding_energies_by_PSSMs(Seq,input_strs,model,peptide_len)

[peptides]=get_peptide_matrix_from_sequence(Seq,peptide_len);
[PSSMs]=get_PSSM_for_str(input_strs,model);
[Energies]=get_gral_HLA_peptides_energies_by_PSSM(peptides,PSSMs);


