function [Modified_Candidates,num]=transfer_info_from_trained_models(str,exp_binding_energies); 

if nargin < 3
   IsBil = 1; 
end

num_pos_candidates = length(exp_binding_energies);
max_str = 175;

dist = zeros(1,num_pos_candidates);

for i=1:num_pos_candidates
    dist(i) = length(find(abs(str(1:max_str)-exp_binding_energies(i).str(1:max_str))));
end

[minm,index]=min(dist);

num = index;

len = exp_binding_energies(num).len;
Modified_Candidates.length = len;


Modified_Candidates.h=exp_binding_energies(num).h;   
Modified_Candidates.pos=exp_binding_energies(num).pos;


Modified_Candidates.thr_pos=exp_binding_energies(num).thr_pos;    
Modified_Candidates.thr_percentages=exp_binding_energies(num).thr_percentages;
Modified_Candidates.set=exp_binding_energies(num).set;

for i=1:len 
   pos = Modified_Candidates.pos{i};  
   candidates_str = str(pos);
   if ~isempty(candidates_str)
       Modified_Candidates.candidates{i}=get_aa_indexes(candidates_str);  
   else
       Modified_Candidates.candidates{i}=[];
   end

   pos = Modified_Candidates.thr_pos{i};
   candidates_str = str(pos);
   if ~isempty(candidates_str)
       Modified_Candidates.thr_candidates{i}=get_aa_indexes(candidates_str);  
   else
       Modified_Candidates.thr_candidates{i}=[];
   end
end


num_cross = length(exp_binding_energies(1).set);
num_h_iter = length(exp_binding_energies(1).set(1).h_iter);
num_iter = length(exp_binding_energies(1).set(1).h_iter(1).bil);

B_type={'gl_B_lin_reg','gl_B_linp_reg'};
Pair_mtrx = {'miy','bet'};
num_types = length(B_type);
num_mtrx = length(Pair_mtrx);

for i=1:num_cross
   for j=1:num_h_iter
      for k=1:num_iter
         for mm=1:num_mtrx
            for tt=1:num_types
               MM =  Pair_mtrx{mm};
               TT = B_type{tt};
               temp=[TT '_' MM];
               eval(['Modified_Candidates.set(i).h_iter(j).bil(k).' temp '= exp_binding_energies(1).set(i).h_iter(j).bil(k).' temp ';']);
            end
         end
      end
   end
end
