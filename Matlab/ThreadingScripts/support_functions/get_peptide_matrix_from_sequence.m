function [peptide_matrix]=get_peptide_matrix_from_sequence(seq_str,peptide_len);

if nargin < 2
   peptide_len = 9; 
end

seq_len=length(seq_str);
peptide_matrix = zeros(peptide_len,seq_len);

for i=1:peptide_len
   j=i-1;
   peptide_matrix(i,1:end-j)=seq_str(i:end);
end

peptide_matrix = (peptide_matrix(:,1:end-peptide_len+1))';
peptide_matrix = char(peptide_matrix);