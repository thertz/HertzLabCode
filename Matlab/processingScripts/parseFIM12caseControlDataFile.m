
function [FIM12struct] = parseFIM12caseControlDataFile(rootPath,fileName)

% case-control dataset file sent by Andrew Dunning

if(~exist('fileName','var'))
	fileName = '/ArrayData/Influenza/FIM12/docs/FIM12_mircoarray_samples_12Dec13.txt';
end


%sub_id	case_control	trta	age_gt_med	SEX_CL


fd = fopen([rootPath,fileName],'r');
[C] = textscan(fd,'%s %s %s %d %s','HeaderLines',1,'Delimiter','\t');
fclose(fd);


for i=1:length(C{1})
	FIM12struct(i).ptid      = C{1}{i}(5:end);
	FIM12struct(i).site      = C{1}{i}(1:3);
	FIM12struct(i).slideID   = [FIM12struct(i).ptid,'_',FIM12struct(i).site];
	FIM12struct(i).gender    = C{5}{i};
	FIM12struct(i).treatment = C{3}{i};
	if(strcmp(FIM12struct(i).treatment,'Fluzone'))
		FIM12struct(i).label = 0;
	else
		FIM12struct(i).label = 1;
	end	
	FIM12struct(i).case_control = C{2}{i};
	if(strcmp(FIM12struct(i).case_control, 'Case'))
		FIM12struct(i).H3N2pos = 1;
	else
		FIM12struct(i).H3N2pos = 0;
	end
end
