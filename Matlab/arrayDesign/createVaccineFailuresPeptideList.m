clear all;
close all;

rootPath = '/Users/thertz/';
loadPath = ['Work/data/AntigenMicroarrays/VaccineFailures/'];

fastaFilename  = [loadPath,'VaccineFailuresInfluenzaStrains.fasta'];
%code is ready, size of peptide and overlap are what determine what we get to see... TBD...
pepLength    = 20;
pepOverlap   = 15;
indexStrains = [1 2];

insertionPositions = { [] [] [] [] [] [] [146] };

[peptideStruct] = createPeptideSetFromFasta(fastaFilename,pepLength,pepOverlap,insertionPositions)



