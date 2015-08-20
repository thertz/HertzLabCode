function [Energies]=get_eng_symmetric_specific_model(exp_binding_energies,B_Glo,B_Spe,peptides,input_matrix);

types = {'gl_linp'};

num_aa = 20;
dim=num_aa*(num_aa+1)/2;   
index_aa=cumsum(1:num_aa);
index_aa = [0 index_aa(1:end-1)];


temp_linp=transpose(B_Glo(1:dim,1)).*input_matrix;

[num_peptides,len] = size(peptides);

num = zeros(1,len);

for j=1:len
    num(j) = length(exp_binding_energies.candidates{j});      
end

num_fin = cumsum(num);
num_ini =[1 num_fin(1:end-1)+1];
num_cand_coeff =num_fin(end);
dim12 = num_cand_coeff.*num_peptides;    


ref = 1:num_peptides;

bi =ones(num_peptides,1)*[exp_binding_energies.h{:}];
aa_indxs =zeros(num_peptides,num_cand_coeff);
for j=1:len       
   temp_cell = exp_binding_energies.candidates{j};
   temp_length=length(temp_cell);
   if temp_length
      ind1 = peptides(:,j)';
      for kk=1:temp_length
         ind2 = temp_cell(kk); 
         inter = num_ini(j)+kk-1;
         aa_indxs(:,inter) = transpose(index_aa(max(ind1,ind2))+min(ind1,ind2));
      end
   end
end

matrix = reshape(temp_linp(reshape(aa_indxs,num_peptides,num_cand_coeff)).*reshape(bi,num_peptides,num_cand_coeff),num_peptides,num_cand_coeff); 
Energies = [matrix ones(num_peptides,1)]*B_Spe + B_Glo(end);   
                 
