% Translate sequence matrix to numerical vectors, according to aa2vec.
% Note: Peptides with X (undeterminced) residues are ignored.
% aa2Vec is an optional parameter, if not specified, uses the VB aa2Vec
% struct...
function [M, validInd] = SequenceToVector(seq, aa2vec)

if(~exist('aa2Vec','var'))
   load('VB.mat');
end

[nSeqs, seqlen] = size(seq);
preLength = size(aa2vec.M, 1);

aaInd = zeros(size(seq));
for i=1:length(aa2vec.aa),
   ind = find(seq==aa2vec.aa{i});
   aaInd(ind) = i;
end

[indX, d] = find(~aaInd);
validInd = setdiff(1:nSeqs, indX);
validSeq = aaInd(validInd, :)';
M = reshape(aa2vec.M(:,validSeq(:)), [seqlen*preLength, size(validSeq,2)])';

validSeq = validSeq';
