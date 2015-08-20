function [indxs_Mtrx1_on_Mtrx2,indxs_Mtrx2_on_Mtrx1]=find_matches_on_str_matrices(Mtrx1,Mtrx2);
%[indxs_Mtrx1_on_Mtrx2,indxs_Mtrx2_on_Mtrx1]=find_matches_on_str_matrices(Mtrx1,Mtrx2);

%max_matrix_size = 30000000;
%max_matrix_size = 5000000;
max_matrix_size = 3000000;


[num1,len1]=size(Mtrx1);
[num2,len2]=size(Mtrx2);

if len1 ~= len2
    error('Str-Mtrxs 1 and 2 should have the same string len');
else
    len=len1;clear len1 len2;
end

total_dim = num1*num2*len;

if  total_dim < max_matrix_size
    compare=sum(abs(repmat(Mtrx1,num2,1)-reshape(repmat((Mtrx2(:))',num1,1),num1*num2,len)),2);
    compare = reshape(compare,num1,num2);
    [indxs_Mtrx1_on_Mtrx2,indxs_Mtrx2_on_Mtrx1]=find(compare==0);
else
    num_blocks = ceil(total_dim/max_matrix_size);
    swap=0;
    if num1 < num2
       swap=1;
       Mtrx2b=Mtrx2;
       Mtrx2=Mtrx1;Mtrx1=Mtrx2b;clear Mtrx2b;
       num2b=num2;
       num2=num1;num1=num2b;clear num2b;
    end

    block_num1 = floor(num1/num_blocks)*ones(1,num_blocks);
    num_left_out = num1-(block_num1(1).*(num_blocks));
    block_num1(1:num_left_out)=block_num1(1:num_left_out)+1;
    block_num1_end = cumsum(block_num1);
    block_num1_ini = [1 block_num1_end(1:end-1)+1];
    indxs_Mtrx1_on_Mtrx2=[];indxs_Mtrx2_on_Mtrx1=[];
    for i=1:num_blocks
       inter =  [block_num1_ini(i):block_num1_end(i)];
       compare=sum(abs(repmat(Mtrx1(inter,:),num2,1)-reshape(repmat((Mtrx2(:))',block_num1(i),1),block_num1(i)*num2,len)),2);
       compare = reshape(compare,block_num1(i),num2);
       [block_indxs_Mtrx1_on_Mtrx2,block_indxs_Mtrx2_on_Mtrx1]=find(compare==0);
       indxs_Mtrx1_on_Mtrx2 = [indxs_Mtrx1_on_Mtrx2 ; block_indxs_Mtrx1_on_Mtrx2+block_num1_ini(i)-1];
       indxs_Mtrx2_on_Mtrx1 = [indxs_Mtrx2_on_Mtrx1 ; block_indxs_Mtrx2_on_Mtrx1];
    end
    
    if swap
        indxs_Mtrx1_on_Mtrx2b = indxs_Mtrx1_on_Mtrx2;
        indxs_Mtrx1_on_Mtrx2 = indxs_Mtrx2_on_Mtrx1;
        indxs_Mtrx2_on_Mtrx1 = indxs_Mtrx1_on_Mtrx2b;
    end
end
