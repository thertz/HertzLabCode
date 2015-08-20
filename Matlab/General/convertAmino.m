function tbl = convertAmino(amino)
aminoN = amino;
for i=1:size(amino,2)
    aminoN{i}{1} = i;
    for l=1:length(amino{i}{4})
        aminoN{i}{4}{l} = convertNucleotide(amino{i}{4}{l});
    end
end

tbl = zeros(4,4,4);
assigned = zeros(1,20);
for i = 1:4
    for j = 1:4
        for k = 1:4
            for l=1:20
                for m=1:length(aminoN{l}{4})
                    if (sum(abs([i j k] - aminoN{l}{4}{m})) == 0)
                        if (tbl(i,j,k) ~= 0)
                            error('update');
                        end
                        tbl(i,j,k) = l;
                        assigned(l) = 1;
                    end
                end
            end
        end
    end
end

if (sum(assigned) ~= 20)
    error('assigned');
end
characteristics = {'hydrophobic' 'positive' 'negative' 'polar' 'small' 'tiny' 'aromatic' 'aliphatic' 'cyclic' 'buried'}; 
%         H P N P C S T A L
%1  A Ala 1 0 0 0 0 1 1 0 0  
%2  C Cys 1 0 0 0 0 1 0 0 0 
%3  D Asp 0 0 1 1 1 1 0 0 0 
%4  E Glu 0 0 1 1 1 0 0 0 0  
%5  F Phe 1 0 0 0 0 0 0 1 0  
%6  G Gly 1 0 0 0 0 1 1 0 0 
%7  H His 1 1 0 1 1 0 0 1 0 
%8  I Ile 1 0 0 0 0 0 0 0 1  
%9  K Lys 1 1 0 1 1 0 0 0 0  
%10 L Leu 1 0 0 0 0 0 0 0 1  
%11 M Met 1 0 0 0 0 0 0 0 0  
%12 N Asn 0 0 0 1 0 1 0 0 0  
%13 P Pro 0 0 0 0 0 1 0 0 0 
%14 Q Gln 0 0 0 1 0 0 0 0 0  
%15 R Arg 0 1 0 1 1 0 0 0 0  
%16 S Ser 0 0 0 1 0 1 1 0 0 
%17 T Thr 1 0 0 1 0 1 0 0 0  
%18 V Val 1 0 0 0 0 1 0 0 1  
%19 W Trp 1 0 0 1 0 0 0 1 0  
%20 Y Tyr 1 0 0 1 0 0 0 1 0  

