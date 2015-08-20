% returns the blosum62 matrix without any amino acid reshuffling.

function Blosum62_Matrix = get_Blosum62_Matrix();

num_aa = 20;
Blosum62_Matrix = zeros(num_aa,num_aa);
i=0;

i=i+1;Blosum62_Matrix(i,:)=[+9 -1 -1 -3 +0 -3 -3 -3 -4 -3 -3 -3 -3 -1 -1 -1 -1 -2 -2 -2];
i=i+1;Blosum62_Matrix(i,:)=[+0 +4 +1 -1 +1 +0 +1 +0 +0 +0 -1 -1 +0 -1 -2 -2 -2 -2 -2 -3];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +4 +1 -1 +1 +0 +1 +0 -0 -0 -1 -0 -1 -2 -2 -2 -2 -2 -3];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +7 -1 -2 -1 -1 -1 -1 -2 -2 -1 -2 -3 -3 -2 -4 -3 -4];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +4 -0 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -2 -3];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +6 -2 -1 -2 -2 -2 -2 -2 -3 -4 -4 -0 -3 -3 -2];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +6 +1 -0 -0 -1 +0 -0 -2 -3 -3 -3 -3 -2 -4];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +6 +2 +0 -1 -2 -1 -3 -3 -4 -3 -3 -3 -4];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +5 +2 +0 +0 +1 -2 -3 -3 -3 -3 -2 -3];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +5 +0 +1 +1 -0 -3 -2 -2 -3 -1 -2];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +8 -0 -1 -2 -3 -3 -2 -1 +2 -2];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +5 +2 -1 -3 -2 -3 -3 -2 -3];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +5 -1 -3 -2 -3 -3 -2 -3];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +5 +1 +2 -2 -0 -1 -1];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +4 +2 +1 -0 -1 -3];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +4 +3 -0 -1 -2];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +4 -1 -1 -3];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +6 +3 +1];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +7 +2];
i=i+1;Blosum62_Matrix(i,:)=[+0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 +0 11];

for i=2:num_aa
    for j=i-1:-1:1
        Blosum62_Matrix(i,j) = Blosum62_Matrix(j,i);
    end
end

% reshuffle = [2 16 17 13 1 6 12 3 4 14 7 15 9 11 8 10 18 5 20 19];

% [x,indexes]=sort(reshuffle);

% temp = zeros(num_aa,num_aa);

% for i=1:num_aa
%    for j=1:num_aa
%        temp(i,j) = Blosum62_Matrix(indexes(i),indexes(j));
%    end
% end

% Blosum62_Matrix = temp;