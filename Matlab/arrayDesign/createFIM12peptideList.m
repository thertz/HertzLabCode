gclear all;
close all;

rootPath = '/Volumes/Data/';
loadPath = ['Work/data/AntigenMicroarrays/FIM12/'];

fastaFilename  = [loadPath,'FIM12_Y2_strainsForArrays.fasta'];
%code is ready, size of peptide and overlap are what determine what we get to see... TBD...
pepLength    = 20;
pepOverlap   = 15;
indexStrains = [1 2 3];

insertionPositions = { [] [] [] [] [] [] [] [] [] [] [] };

[peptideStruct] = createPeptideSetFromFasta(fastaFilename,pepLength,pepOverlap,insertionPositions,rootPath)



