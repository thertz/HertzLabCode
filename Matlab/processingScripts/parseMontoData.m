
function [ptidStruct] = parseMontoData(rootPath,fileName)

if(~exist('fileName','var'))
  fileName = 'ArrayData/Influenza/Monto/Tomer_UM_sera0708_HAI_titers.txt';
end

fd = fopen([rootPath,fileName],'r');

%StudyID	AGE	ELIG_ILLNESS	TREATMENT	VAX_DATE	VAXCODE2	S9_date	S9_ACCN	S10_date	S10_ACCN	S11_date	S11_ACCN	AH3_WIv_s9t	AH3_WIv_s10t	AH3_WIv_s11t	AH3_UGc_s9t	AH3_UGc_s10t	AH3_UGc_s11t
[C] = textscan(fd,'%s %d %d %s %s %d %s %s %s %s %s %s %d %d %d %d %d %d %d','Delimiter','\t','HeaderLines',1);
fclose(fd);

for i=1:length(C{1})

	ptidStruct(i).ptid                = C{1}{i};
	ptidStruct(i).age                 = C{2}(i);
	ptidStruct(i).treatmentOrig       = C{4}{i};
	ptidStruct(i).H3_Wis_HI_preTiter  = C{13}(i); 
	ptidStruct(i).H3_Wis_HI_postTiter = C{14}(i); 
	ptidStruct(i).H3_Wis_HI_ltTiter   = C{15}(i); 
	ptidStruct(i).H3_Ug_HI_preTiter   = C{16}(i); 
	ptidStruct(i).H3_Ug_HI_postTiter  = C{17}(i); 
	ptidStruct(i).H3_Ug_HI_ltTiter    = C{18}(i); 
end