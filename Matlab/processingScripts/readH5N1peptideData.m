
% read mapping of 96 wells of original peptide set of 18mers from St. Jude
% from June 2010.

function [pepData] = readH5N1peptideData(rootPath)

if(~exist('rootPath','var'))
  rootPath = '/Users/thertz/';
end

fileName = 'Documents/Projects/AntigenMicroarrays/Influenza/H5N1/H5_Avian_HA_peptides_L20_O10_18merMapping_June_28_10.txt'; 
 
fd = fopen([rootPath,fileName],'r');

[C] = textscan(fd,'%s %d %d %s %s %d','Delimiter','\t','HeaderLines',1);
fclose(fd);
 
for i=1:length(C{1})
  pepData(i).strain   = C{1}{i};
  pepData(i).begInd   = C{2}(i)+2;
  pepData(i).endInd   = C{3}(i);
  pepData(i).sequence = C{5}{i};
  pepData(i).tubeNo   = C{6}(i);
end
