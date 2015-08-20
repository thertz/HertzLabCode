
function [FIM12struct] = parseFIM12dataFile(rootPath,fileName)

if(~exist('fileName','var'))
	fileName = '/ArrayData/Influenza/FIM12/docs/fim12.csv';
end

%SUB_ID	SEX_CL	frail_conds	comorb_conds	stsp_indx	mitnit	trta	trtp	FullAS	PPAS	vacdate	dobdate	flu_vac_prec	sub_id_yr1	
%termdate	Startdate	PD_ILI	LabConf	subtyp_lineg


fd = fopen([rootPath,fileName],'r');
[C] = textscan(fd,'%s %s %d %d %f %f %s %s %d %d %s %s %d %s %s %s %d %d %s\n','HeaderLines',1,'Delimiter',',');
fclose(fd);


for i=1:length(C{1})
	FIM12struct(i).ptid      = C{1}{i}(5:end);
	FIM12struct(i).site      = C{1}{i}(1:3);
	FIM12struct(i).slideID   = [FIM12struct(i).ptid,'_',FIM12struct(i).site];
	FIM12struct(i).gender    = C{2}{i};
	FIM12struct(i).treatment = C{7}{i};
	if(strcmp(FIM12struct(i).treatment,'Fluzone'))
		FIM12struct(i).label = 0;
	else
		FIM12struct(i).label = 1;
	end	
	FIM12struct(i).PD_ILI    = C{17}(i);
	FIM12struct(i).LabConf   = C{18}(i);
	FIM12struct(i).subtype   = C{19}(i);
	if(strcmp(FIM12struct(i).subtype,'H3N2'))
		FIM12struct(i).H3N2pos = 1;
	else
		FIM12struct(i).H3N2pos = 0;
	end
end
