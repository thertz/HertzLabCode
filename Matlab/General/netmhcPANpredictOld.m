
%run NETMHCpan 2.0 standalone on a fasta file using a list of alleles
function [] = netmhcPANpredict(seqFile,alleleList,resName,resPath,pepLength,rootPath)

if(~exist('rootPath','var'))
  rootPath = '/Volumes/E/';
end

if(~exist('pepLength','var'))
  pepLength = 9;
end

allAlleles = [];
for i=1:length(alleleList)

  %use required allele format...
  currAllele = strrep(alleleList{i},'_','');
  currAllele = strrep(currAllele,'*','');
  currAllele = ['HLA-',currAllele];
  
  if(i==1)
    allAlleles = [currAllele];
  else
    allAlleles = [allAlleles,',',currAllele];
  end
end

resFile = [rootPath,resPath,resName,'_NETMHCpanPredictions.txt'];
if(exist(resFile,'file'))
  system([rootPath,'Work/packages/netMHCpan-2.0/netMHCpan -a ',allAlleles,' -l ',num2str(pepLength),' ',seqFile,' >> ' resFile]);
else
  system([rootPath,'Work/packages/netMHCpan-2.0/netMHCpan -a ',allAlleles,' -l ',num2str(pepLength),' ',seqFile,' > ' resFile]);
end
