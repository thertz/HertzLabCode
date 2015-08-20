
function [sampleDateStruct] = parseFIM12SamplesByDates(rootPath,fileName)

	
if(~exist('fileName','var'))
	fileName = 'ArrayData/Influenza/FIM12/docs/FIM12samplesByDates.txt';
end

%sample ID	slide	Date

fd = fopen([rootPath,fileName],'r');
[C] = textscan(fd,'%s %d %s','HeaderLines',1,'Delimiter','\t');
fclose(fd);


for i=1:length(C{1})
	sampleDateStruct(i).sampleID = C{1}{i};
	sampleDateStruct(i).slide    = C{2}(2);
	sampleDateStruct(i).date     = C{3}{i};
end

