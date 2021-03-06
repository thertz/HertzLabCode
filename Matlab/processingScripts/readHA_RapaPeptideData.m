
% read mapping of 96 wells of original peptide set of 18mers from St. Jude
% from Jan 2012.

function [pepData] = readHA_RapaPeptideData(rootPath)

if(~exist('rootPath','var'))
  rootPath = '/Volumes/viddshared/HertzLab/';
end

% headers is:
% Strain	beg Ind	end Ind	sequence Synthesized Sequence	Peptide Number	Printing Number


fileName = ['ArrayData',filesep,'Influenza',filesep,'Rapa',filesep,'docs',filesep,'McGargillStrainsFinal_peptides_L20_O15.txt'];

fd = fopen([rootPath,fileName],'r');

[C] = textscan(fd,'%s %d %d %s %s %d %d','Delimiter','\t','HeaderLines',1);
fclose(fd);
 
for i=1:length(C{1})
  pepData(i).strain   = C{1}{i};
  pepData(i).begInd   = C{2}(i);
  pepData(i).endInd   = C{3}(i);
  pepData(i).sequence = C{4}{i};
  
  %sanity check synthesized seq is identical to ordered seq:
  if(~strcmp(pepData(i).sequence,C{5}{i}))
    error(sprintf('syntesized sequence is not identical to ordered sequence:\n%s\n%s\n',C{4}{i},C{5}{i}));
  end
  pepData(i).printedInd = C{7}(i);
end
