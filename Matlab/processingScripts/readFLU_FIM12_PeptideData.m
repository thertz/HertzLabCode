
% read mapping of peptides for Dengue project based on Nicaragua strains - Eva Harris
% from July 2013. Ryan constructed this mapping and verified that peptide ordering did NOT change. Control sequences
% were not synthesized due to error on their end.

function [pepData] = readFLU_FIM12_PeptideData(rootPath)

% headers is:
%Printed Array Number	SEQUENCE	Plate Row	Plate Column	Strain	beg Ind	end Ind

fileName = ['ArrayData',filesep,'Influenza',filesep,'FIM12',filesep,'docs',filesep,'FIM12_Strains_PeptidesPrinted_May2014.txt'];

fd = fopen([rootPath,fileName],'r');

[C] = textscan(fd,'%s %s %s %s %d %d','Delimiter','\t','HeaderLines',1);
fclose(fd);
 
for i=1:length(C{1})
  pepData(i).name     = C{1}{i};
  pepData(i).sequence = C{2}{i};
  pepData(i).modifiedSequence = C{3}{i};
  pepData(i).strain   = C{4}{i};
  pepData(i).begInd   = C{5}(i);
  pepData(i).endInd   = C{6}(i);
end
