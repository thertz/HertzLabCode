
% read mapping of peptides for Dengue project based on Nicaragua strains - Eva Harris
% from July 2013. Ryan constructed this mapping and verified that peptide ordering did NOT change. Control sequences
% were not synthesized due to error on their end.

function [pepData] = readDengueNicaragua_PeptideData(rootPath)

% headers is:
%Printed Array Number	SEQUENCE	Plate Row	Plate Column	Strain	beg Ind	end Ind

fileName = ['ArrayData',filesep,'Dengue',filesep,'docs',filesep,'DenguePeptides_Env_NS1_prM_PrintedArrayList.txt'];

fd = fopen([rootPath,fileName],'r');

[C] = textscan(fd,'%s %s %c %d %s %d %d','Delimiter','\t','HeaderLines',1);
fclose(fd);
 
for i=1:length(C{1})
  pepData(i).name     = C{1}{i};
  pepData(i).sequence = C{2}{i};
  pepData(i).strain   = C{5}{i};
  pepData(i).begInd   = C{6}(i);
  pepData(i).endInd   = C{7}(i);
end
