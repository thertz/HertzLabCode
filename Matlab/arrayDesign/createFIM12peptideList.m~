clear all;
close all;

rootPath = '/Volumes/Data/';
loadPath = ['Work/data/AntigenMicroarrays/H7N9/'];

fastaFilename  = [loadPath,'H7N9_strainsForArrays.fasta'];
%code is ready, size of peptide and overlap are what determine what we get to see... TBD...
pepLength    = 20;
pepOverlap   = 10;
indexStrains = [1 2 3];

insertionPositions = { [] [] [] [] [] [] [] [] [] [] [] };

[peptideStruct] = createPeptideSetFromFasta(fastaFilename,pepLength,pepOverlap,insertionPositions,rootPath)



