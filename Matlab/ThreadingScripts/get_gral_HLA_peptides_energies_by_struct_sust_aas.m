function [Energies]=get_gral_HLA_peptides_energies_by_struct_sust_aas(peptides,HLA_Strings,model,pos,aas,cross_set,iter_h,iter);
%[Energies]=get_gral_HLA_peptides_energies_by_struct_sust_aas(peptides,HLA_Strings,model,aas,pos,cross_set,iter_h,iter);
if nargin < 8
    iter=length(model(1).set(1).h_iter.bil);
    if nargin < 7
        iter_h=1;
        if nargin < 6
            cross_set=1;
        end
    end
end

peptide_len = size(peptides,2);

num_aas = length(aas);
num_pos = length(pos);

if num_pos ~= num_aas
    error('Number of position to substitute in pos by the number of aminoacids in aa do not match');
end

if ischar(peptides)
   peptides = get_aa_indexes(peptides);    
end


num_HLAs = length(HLA_Strings);

num_peptides = size(peptides,1);

Energies=zeros(num_peptides,num_HLAs);

num_aa = 20;
dim1=num_aa*(num_aa+1)/2 + 1; 
dim2 = num_aa*num_aa + 1;

input_matrix = model(1).matrix;
for i=1:num_HLAs
   HLA_str = HLA_Strings(i).str;   
   [Modified_Candidates,num]=transfer_info_from_trained_uni_models(HLA_str,model); 
   HLA_str(pos)=aas;
   [Modified_Candidates]=transfer_info_from_trained_uni_models(HLA_str,model,num);
%   fprintf(['HLA ' HLA_name{i} ' has as its closest model. The model for HLA : ' exp_binding_energies(num).name '\n']);  
   [B_Struct] = get_B_struct_from_uni_exp_binding_structs(Modified_Candidates,cross_set,iter_h,iter);           
   vec_len = length(B_Struct.Glo.miy.SS_all.gl_linp);
   
   
   B_Glo = B_Struct.Glo.miy.SS_all.gl_linp;
   B_Spe = B_Struct.Spe(1).miy.SS_all.gl_linp;

   if vec_len == dim2
       [Energies(:,i)]=-get_eng_both_ways_specific_model(Modified_Candidates,B_Glo,B_Spe,peptides,input_matrix);
   else
       if vec_len == dim1
           [Energies(:,i)]=-get_eng_symmetric_specific_model(Modified_Candidates,B_Glo,B_Spe,peptides,input_matrix);
       else
           error(['Can support a model with ' num2str(vec_len) ' weights ; ']); 
       end
   end
end
