function nmbr = convert(nucleotide)
alphabet = 'ACTG';
empty = '- ';
for i=1:length(nucleotide)
    a = find(upper(nucleotide(i)) == alphabet);
    e = find(nucleotide(i) == empty);
    %assert(~isempty(a) | ~isempty(e));
    if (isempty(e))
        nmbr(i) = a;
    else
        nmbr(i) = 0;
    end
end
