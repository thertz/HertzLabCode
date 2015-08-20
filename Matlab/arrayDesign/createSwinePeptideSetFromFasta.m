clear all;
close all;

rootPath = '/Volumes/Data/';
loadPath = ['Work/data/AntigenMicroarrays/SwineFlu/'];

fastaFilename  = [loadPath,'H1_H3_SwineStrains_CrystalLoving.fasta'];
%Code is ready, size of peptide and overlap are what determine what we get to see... TBD...
pepLength    = 20;
pepOverlap   = 15;
indexStrains = [1:2];

insertionPositions = { [] [] []};

[peptideStruct] = createPeptideSetFromFasta(fastaFilename,pepLength,pepOverlap,insertionPositions,rootPath)

fd = fopen([rootPath,loadPath,'H1_H3_SwineStrains_peptideList_20mer_15overlap'],'w');

fprintf(fd,'Strain\tbegInd\tendInd\tseqquence\n');
for i=1:length(peptideStruct)
  fprintf(fd,'%s\t%d\t%d\t%s\t\n',peptideStruct(i).Header,peptideStruct(i).begInd,peptideStruct(i).endInd,peptideStruct(i).sequence);
end
fclose(fd)
