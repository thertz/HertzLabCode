function [B_Struct] = get_B_struct_from_uni_exp_binding_structs(exp_binding_energies,cross_set,h_iter_num,iter_num);

Pair_mtrx = {'miy','bet'};
specific_regression_type = {'linp_reg_linp_reg'};
global_regression_type = {'linp_reg'};
B_type={'gl_linp'};
%specific_regression_type = {'linp_reg_lin_reg','linp_reg_linp_reg'};
%global_regression_type = {'lin_reg','linp_reg'};
%B_type={'gl_lin','gl_linp'};

num_types = length(B_type);
num_mtrx = length(Pair_mtrx);

if nargin>2
   for mm=1:num_mtrx
      MM = Pair_mtrx{mm};
      num_HLA_types = length(exp_binding_energies);
      for tt = 1:num_types
         TT = B_type{tt};TT_Spe = specific_regression_type{tt};TT_Glo = global_regression_type{tt}; 
         if iter_num > 1
             eval(['B_Struct.Glo.' MM '.SS_all.' TT ' = exp_binding_energies(1).set(cross_set).h_iter(h_iter_num).bil(iter_num-1).gl_B_' TT_Glo '_' MM ';']);
         else
             if tt == 1
                 eval(['B_Struct.Glo.' MM '.SS_all.' TT ' = ones(size(exp_binding_energies(1).set(cross_set).h_iter(h_iter_num).bil(iter_num).gl_B_' TT_Glo '_' MM '));']);
             else
                 eval(['B_Struct.Glo.' MM '.SS_all.' TT ' = ones(size(exp_binding_energies(1).set(cross_set).h_iter(h_iter_num).bil(iter_num).gl_B_' TT_Glo '_' MM '));']);
                 eval(['B_Struct.Glo.' MM '.SS_all.' TT '(end)=0;']);
             end
         end
         for k=1:num_HLA_types
            eval(['B_Struct.Spe(k).' MM '.SS_all.'  TT ' = exp_binding_energies(1).set(cross_set).h_iter(h_iter_num).bil(iter_num).B_' TT_Spe '_' MM ';']); 
         end
      end
   end
else
   [num_cross]=length(exp_binding_energies(1).set);
   [num_iter_h]=length(exp_binding_energies(1).set(1).h_iter);
   [num_iter]=length(exp_binding_energies(1).set(1).h_iter(1).bil);
   for mm=1:num_mtrx
      MM = Pair_mtrx{mm};
      eval(['load ' str MM ';']);
      
      num_HLA_types = length(exp_binding_energies)-1;
   
      for tt = 1:num_types
         for cross_set=1:num_cross
            for h_iter_num=1:num_iter_h
               for iter_num=1:num_iter
                  TT = B_type{tt};TT_Spe = specific_regression_type{tt};TT_Glo = global_regression_type{tt}; 
                  eval(['B_Struct(cross_set,iter_num).Glo.' MM '.SS_all.' TT ' = exp_binding_energies(1).set(cross_set).h_iter(h_iter_num).bil(iter_num).gl_B_' TT_Glo '_' MM ';']);
                  for k=1:num_HLA_types
                     eval(['B_Struct(cross_set,iter_num).Spe(k).' MM '.SS_all.' TT ' = exp_binding_energies(1).set(cross_set).h_iter(h_iter_num).bil(iter_num).B_' TT_Spe '_' MM ';']); 
                  end
               end
            end
         end
      end
   end
end

