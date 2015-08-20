function [tbl,letter,nmbr,labels,amino] = AminoTables()
[amino,letter,nmbr,labels] = buildAmino;
tbl = convertAmino(amino);

