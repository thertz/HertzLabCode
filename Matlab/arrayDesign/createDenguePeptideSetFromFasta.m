clear all;
close all;

rootPath = '/Users/thertz/';
loadPath = ['Work/data/AntigenMicroarrays/Dengue/'];

fastaFilename  = [loadPath,'NicaraguaDengueStrains_Env_NS1_prM.fasta'];
%Code is ready, size of peptide and overlap are what determine what we get to see... TBD...
pepLength    = 20;
pepOverlap   = 15;
indexStrains = [1:9];

insertionPositions = { [] [] [] [] [] [] [] [] [] };

[peptideStruct] = createPeptideSetFromFasta(fastaFilename,pepLength,pepOverlap,insertionPositions)

fd = fopen([rootPath,loadPath,'NicaraguaDengueStrains_Env_NS1_prM_peptideList_20mer_15overlap'],'w');

fprintf(fd,'Strain\tbegInd\tendInd\tseqquence');
for i=1:length(peptideStruct)
  fprintf(fd,'%s\t%d\t%d\t%s\t\n',peptideStruct(i).Header,peptideStruct(i).begInd,peptideStruct(i).endInd,peptideStruct(i).sequence);
end
fclose(fd)
