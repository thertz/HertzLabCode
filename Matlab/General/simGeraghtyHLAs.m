
clear  all;
close all;
%loadPath = '/Volumes/E/Work/Results/InfluenzaStrains/'
%fileName = 'ThomasLab-HLATypingReport-073010.txt';

loadPath = '~/Work/Results/';
fileName = 'ThomasLab-HLATypingReportCombined.txt';


%PROTOCOL	LABID	ASSAYID	SPECID	LIMSID	PTID	VISITNO	DRAWDT	ASSAYRUN	HLALOCUS	GENOTYPE	GTAPPROX	IMGTVER	DIPAMBIG	TESTDT	RELIABLE	REPLACE	MODDT	COMMENTS

fd = fopen([loadPath,fileName],'r');
[C] = textscan(fd,'%d %s %s %s %s %s %d %s %c %s %s %c %s %c %s %c %s %s %s','delimiter','\t','headerLines',1); 
fclose(fd);

ptids = C{6};
HLAs = C{11};

uniqPtids = unique(ptids);

for i=1:length(uniqPtids)
  currInds = strmatch(uniqPtids{i},ptids);

  HLAstruct(i).HLAs    = {};
  HLAstruct(i).rawHLAs = {};
  for j=1:length(currInds)
    HLAstruct(i).ptid = uniqPtids{i};
    currHLAs = decodeGeraghtyHLAs(HLAs(currInds(j)));
    HLAstruct(i).HLAs = [HLAstruct(i).HLAs currHLAs];
    HLAstruct(i).rawHLAs = [HLAstruct(i).rawHLAs repmat(HLAs(currInds(j)),1,length(currHLAs))];
  end
end

fd = fopen([loadPath,fileName(1:end-4),'_decoded.txt'],'w');
fprintf(fd,'%s\t%s\t%s\t%s\n','PTID','HLALOCUS','HLA','GENOTYPE')
for i=1:length(HLAstruct)
  for j=1:length(HLAstruct(i).HLAs)
    sepInd = findstr('_',HLAstruct(i).HLAs{j});
    fprintf(fd,'%s\t%s\t%s\t%s\n',HLAstruct(i).ptid,['HLA-',HLAstruct(i).HLAs{j}(1:sepInd-1)],strrep(HLAstruct(i).HLAs{j},'_','*'),HLAstruct(i).rawHLAs{j});
  end
end