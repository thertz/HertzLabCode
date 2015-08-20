function [PSSMs]=get_PSSM_for_HLA(HLA_Strings,model,cross_set,iter_h,iter)

if nargin < 5
    iter=length(model(1).set(1).h_iter.bil);
    if nargin < 4
        iter_h=1;
        if nargin < 3
            cross_set=1;
        end
    end
end

max_pos_len = 175;
num_aa=20;
len = model(1).len;

num_model_HLAs = length(model);
model_strs = zeros(num_model_HLAs,max_pos_len);

for i=1:num_model_HLAs
   model_strs(i,:)= model(i).str(1:max_pos_len);   
end

num_input_HLAs = length(HLA_Strings);

input_strs = zeros(num_input_HLAs,max_pos_len);

for i=1:num_input_HLAs
    input_strs(i,:)= HLA_Strings(i).str(1:max_pos_len);   
end
input_strs=char(input_strs);

input_aas = get_aa_indexes(input_strs);

dist = sum(abs(repmat(model_strs,num_input_HLAs ,1)-reshape(permute(reshape(repmat(input_strs,1,num_model_HLAs),num_input_HLAs,max_pos_len,num_model_HLAs),[3,1,2]),num_model_HLAs*num_input_HLAs,max_pos_len)),2);
dist =reshape(dist,num_model_HLAs,num_input_HLAs)';

[x,num] = min(dist,[],2);

[B_Struct] = get_B_struct_from_uni_exp_binding_structs(model(1),cross_set,iter_h,iter);

PSSMs = cell(1,num_input_HLAs);

dim = num_aa*num_aa;
GL_matrix = (B_Struct.Glo.miy.SS_all.gl_linp(1:end-1))';
GL_matrix = GL_matrix.*model(1).matrix;

if length(GL_matrix)== dim
    GL_matrix = reshape(GL_matrix,num_aa,num_aa);
else
    GL_matrix = inverse_compact_matrices(GL_matrix);
end
GL_matrix = GL_matrix';

W = B_Struct.Spe(1).miy.SS_all.gl_linp(1:end-1);
W_per_pos = cell(1,len);
so_far=0;
for i=1:len
    num_cand =length(model(1).h{i});
    W_per_pos{i}=W(so_far+1:so_far+num_cand);
    so_far = so_far+num_cand;
end

onesVec = ones(1,num_aa);
for i=1:num_input_HLAs
   PSSMs{i}=zeros(len,num_aa+1);
   for j=1:len
      PSSMs{i}(j,1:num_aa) = sum(GL_matrix(input_aas(i,model(1).pos{j}),:).*((model(num(i)).h{j}'.*W_per_pos{j})*onesVec));
   end
   PSSMs{i}(:,num_aa+1) = B_Struct.Glo.miy.SS_all.gl_linp(end)+B_Struct.Spe(1).miy.SS_all.gl_linp(end);
   PSSMs{i} = PSSMs{i}';
end






