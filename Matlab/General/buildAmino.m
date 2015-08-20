function [amino,letter,nmbr,labels] = buildAmino()
letter = zeros(1,128);
nmbr = zeros(1,20);
amino = cell(1,25);
labels = zeros(10,20);
% wikipeida
% amino{1} = {'A'    'Ala'      'Alanine'          {'GCT' 'GCC' 'GCA'  'GCG'}             [1 0 0 0 0 1 1 0 0]};
% amino{2} = {'C'    'Cys'      'Cysteine'         {'TGT' 'TGC'}                          [1 0 0 0 0 1 0 0 0]};
% amino{3} = {'D'    'Asp'      'Aspartic'         {'GAT' 'GAC'}                          [0 0 1 1 1 1 0 0 0]};
% amino{4} = {'E'    'Glu'      'Glutamic'         {'GAA' 'GAG'}                          [0 0 1 1 1 0 0 0 0]};
% amino{5} = {'F'    'Phe'      'Phenylalanine'    {'TTT' 'TTC'}                          [1 0 0 0 0 0 0 1 0]};
% amino{6} = {'G'    'Gly'      'Glycine'          {'GGT' 'GGC' 'GGA' 'GGG'}              [1 0 0 0 0 1 1 0 0]};
% amino{7} = {'H'    'His'      'Histidine'        {'CAT' 'CAC'}                          [1 1 0 1 1 0 0 1 0]};
% amino{8} = {'I'    'Ile'      'Isoleucine'       {'ATT' 'ATC' 'ATA'}                    [1 0 0 0 0 0 0 0 1]};
% amino{9} = {'K'    'Lys'      'Lysine'          {'AAA' 'AAG'}                           [1 1 0 1 1 0 0 0 0]};
% amino{10} = {'L'   'Leu'      'Leucine'          {'TTG' 'TTA' 'CTT' 'CTC' 'CTA' 'CTG'}  [1 0 0 0 0 0 0 0 1]};
% amino{11} = {'M'   'Met'      'Methionine'       {'ATG'}                                [1 0 0 0 0 0 0 0 0]};
% amino{12} = {'N'   'Asn'      'Asparagine'     {'AAT' 'AAC'}                            [0 0 0 1 0 1 0 0 0]};
% amino{13} = {'P'   'Pro'      'Proline'        {'CCT' 'CCC' 'CCA' 'CCG'}                [0 0 0 0 0 1 0 0 0]};
% amino{14} = {'Q'   'Gln'      'Glutamine'       {'CAA' 'CAG'}                           [0 0 0 1 0 0 0 0 0]};
% amino{15} = {'R'   'Arg'      'Arginine'         {'CGT' 'CGC' 'CGA' 'CGG' 'AGA' 'AGG'}  [0 1 0 1 1 0 0 0 0]};
% amino{16} = {'S'   'Ser'      'Serine'           {'TCT' 'TCC' 'TCA' 'TCG' 'AGT' 'AGC'}  [0 0 0 1 0 1 1 0 0]};
% amino{17} = {'T'   'Thr'      'Threonine'        {'ACT' 'ACC' 'ACA' 'ACG'}              [1 0 0 1 0 1 0 0 0]};
% amino{18} = {'V'   'Val'      'Valine'          {'GTT' 'GTC' 'GTA' 'GTG'}               [1 0 0 0 0 1 0 0 1]};
% amino{19} = {'W'   'Trp'      'Tryptophan'       {'TGG'}                                [1 0 0 1 0 0 0 1 0]};
% amino{20} = {'Y'   'Tyr'      'Tyrosine'         {'TAT' 'TAC'}                          [1 0 0 1 0 0 0 1 0]};


%                                                                                        H*P*N*P*M*S*A*L*C*B
amino{1} = {'A'    'Ala'      'Alanine'          {'GCT' 'GCC' 'GCA'  'GCG'}             [1 0 0 0 0 1 0 1 0 1]};

amino{2} = {'C'    'Cys'      'Cysteine'         {'TGT' 'TGC'}                          [0 0 0 1 1 0 0 0 0 1]};

amino{3} = {'D'    'Asp'      'Aspartic'         {'GAT' 'GAC'}                          [0 0 1 1 1 0 0 0 0 0]};

amino{4} = {'E'    'Glu'      'Glutamic'         {'GAA' 'GAG'}                          [0 0 1 1 0 0 0 0 0 0]};

amino{5} = {'F'    'Phe'      'Phenylalanine'    {'TTT' 'TTC'}                          [1 0 0 0 0 0 1 0 1 1]};

amino{6} = {'G'    'Gly'      'Glycine'          {'GGT' 'GGC' 'GGA' 'GGG'}              [1 0 0 0 0 1 0 1 0 0]};

amino{7} = {'H'    'His'      'Histidine'        {'CAT' 'CAC'}                          [0 1 0 1 0 0 1 0 1 0]};

amino{8} = {'I'    'Ile'      'Isoleucine'       {'ATT' 'ATC' 'ATA'}                    [1 0 0 0 0 0 0 1 0 1]};

amino{9} = {'K'    'Lys'      'Lysine'          {'AAA' 'AAG'}                           [0 1 0 1 0 0 0 0 0 0]};

amino{10} = {'L'   'Leu'      'Leucine'          {'TTG' 'TTA' 'CTT' 'CTC' 'CTA' 'CTG'}  [1 0 0 0 0 0 0 1 0 1]};

amino{11} = {'M'   'Met'      'Methionine'       {'ATG'}                                [1 0 0 0 0 0 0 0 0 1]};

amino{12} = {'N'   'Asn'      'Asparagine'     {'AAT' 'AAC'}                            [0 0 0 1 1 0 0 0 0 0]};

amino{13} = {'P'   'Pro'      'Proline'        {'CCT' 'CCC' 'CCA' 'CCG'}                [1 0 0 0 1 0 0 0 1 0]};

amino{14} = {'Q'   'Gln'      'Glutamine'       {'CAA' 'CAG'}                           [0 0 0 1 0 0 0 0 0 0]};

amino{15} = {'R'   'Arg'      'Arginine'         {'CGT' 'CGC' 'CGA' 'CGG' 'AGA' 'AGG'}  [0 1 0 1 0 0 0 0 0 0]};

amino{16} = {'S'   'Ser'      'Serine'           {'TCT' 'TCC' 'TCA' 'TCG' 'AGT' 'AGC'}  [0 0 0 1 0 1 0 0 0 0]};

amino{17} = {'T'   'Thr'      'Threonine'        {'ACT' 'ACC' 'ACA' 'ACG'}              [0 0 0 1 1 0 0 0 0 0]};

amino{18} = {'V'   'Val'      'Valine'          {'GTT' 'GTC' 'GTA' 'GTG'}               [1 0 0 0 1 0 0 1 0 1]};

amino{19} = {'W'   'Trp'      'Tryptophan'       {'TGG'}                                [1 0 0 0 0 0 1 0 1 1]};

amino{20} = {'Y'   'Tyr'      'Tyrosine'         {'TAT' 'TAC'}                          [1 0 0 0 0 0 1 0 1 0]};

for i =1:20
    c = amino{i}{1};
    letter(double(upper(c))) = i;
    letter(double(lower(c))) = i;
    nmbr(i) = upper(c);
    labels(:,i) = amino{i}{5}';
end

amino{21} = {'B'   'Asx'   'Aspartic/Asparagine'       {'GAT' 'GAC' 'AAT' 'AAC'}};
amino{22} = {'Z'   'Glx'   'Glutamic/Glutamine'       {'GAA' 'GAG' 'CAA' 'CAG'}};
amino{23} = {'\'   'End'       'Terminator'       {'TAA' 'TAG' 'TGA'}};
amino{24} = {'/'   'Beg'       'Start'    {'ATG'}};
amino{25} = {'X'   'Xxx'       'Unknown'          {}};
      