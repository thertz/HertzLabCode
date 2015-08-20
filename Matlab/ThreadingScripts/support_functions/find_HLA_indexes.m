function indexes = find_HLA_indexes(All_HLA_names,HLA_name);

num_targets = length(HLA_name);

uni_length=size(All_HLA_names,2);

indexes = zeros(1,num_targets);

HLA_targets = zeros(num_targets,uni_length);

for i=1:num_targets
    HLA_targets(i,1:length(HLA_name{i})+4)=['HLA_' HLA_name{i}];
end

HLA_targets = char(HLA_targets);
[indxs_Mtrx1_on_Mtrx2,indxs_Mtrx2_on_Mtrx1]=find_matches_on_str_matrices(HLA_targets,All_HLA_names);

indexes(indxs_Mtrx1_on_Mtrx2)=indxs_Mtrx2_on_Mtrx1;

not_found = find(indexes==0);

if ~isempty(not_found)
    warning('Aminoacid Sequences for the following HLAs could not be found: ');
    for i=1:length(not_found)
       fprintf('HLA : %s at index %d on cell array HLA_name\n',HLA_name{not_found(i)},not_found(i)); 
    end
end
