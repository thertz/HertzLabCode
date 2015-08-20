

% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

%
clear all;
close all;

ZachFlag  = 1;
plotFlag  = 0;
printFlag = 0;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'IDRI';

if(ZachFlag)
  rootPath = [filesep,'\fhdata.fhcrc.org',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
else % Tomer path's
  %rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  rootPath = [filesep,'Volumes',filesep,'Data',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
end

addpath([rootPath,'Code',filesep,'General']);
addpath([rootPath,'Code',filesep,'plotScripts']);
  
%multiple dates can be entered, but must have a corresponding prefix
experimentDates = {'05_08_2015','05_13_2015','05_15_2015'};  %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'IDRI_05_08_15_','IDRI_05_13_15_','IDRI_05_15_15_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix



fileSuffix = ''; % deprcated option to have a suffix to the name - we no longer allow this in the current naming format.

%-------------------------------------------------------------------------------%
% Additional parameters - project specific
%-------------------------------------------------------------------------------%
numArrays                = 2; %number of arrays on each slide.
typeFlag                 = 'median'; % uses the median over replicates of an antigen. Can also be 'mean'
colorFlags               = {'635'};
colorTags                = {'IgG'};
switchIDs_and_Names_Flag = 0; % old flag for problem with Scrips GPR files, deprecated
commentFlag              = 1; % slide to sample mapping contain a comment column. Current default, used for backward compatbility
fontSize                 = 20;
%expGroups                = {'<3>','<2>','<1>'};% groups listed in sample to slide mapping
%yLims = [log2(1024) log2(65536)];
%minThreshold             = log2(2048); %minimal threshold for peak responses, used by findPeaks.
yLims        = [0 65000];
minThreshold = 2048; %minimal threshold for peak responses, used by findPeaks.





%-------------------------------------------------------------------------------%
%defining elements outside of loop to add loop data too
slideNames      ={};
ptids           ={};
groupNames      ={};
mPeptResponses  =[];
%-------------------------------------------------------------------------------

for i=1:length(experimentDates)
  
  % name of file in which converted data will be stored for future use
  arrayDataFilename = [savePath,experimentDates{i},filesep,expPrefix{i},'arrayData_',typeFlag,'.mat']; 

  load(arrayDataFilename);
 
  slideNames            = [slideNames arrayData.slideNames];
  ptids                 = [ptids arrayData.ptids];
  groupNames            = [groupNames arrayData.groupNames];
  mPeptResponses        = [mPeptResponses ;[arrayData(1).responseMatrix{1}]];

end
figInd = 1;


