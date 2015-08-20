clear all;
close all;

rootPath = '/Users/thertz/';
loadPath = ['Work/data/AntigenMicroarrays/Rapa/'];

fastaFilename  = [loadPath,'McGargillStrainsFinal.fasta'];
%code is ready, size of peptide and overlap are what determine what we get to see... TBD...
pepLength    = 20;
pepOverlap   = 15;
indexStrains = [1 2 3];

insertionPositions = { [] [] [] };

[peptideStruct] = createPeptideSetFromFasta(fastaFilename,pepLength,pepOverlap,insertionPositions)



