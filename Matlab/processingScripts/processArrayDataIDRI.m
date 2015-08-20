
% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.


clear all;
close all;


ZachFlag = 1;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'IDRI';

if(ZachFlag)
  rootPath = [filesep,filesep,'fhdata.fhcrc.org',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
else % Tomer path's
  %rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  rootPath = [filesep,'Volumes',filesep,'Data',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
end

experimentDates = {'05_13_2015','05_15_2015'};  %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'IDRI_05_13_15_','IDRI_05_15_15_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

fileSuffix = ''; % deprcated option to have a suffix to the name - we no longer allow this in the current naming format.

%-------------------------------------------------------------------------------%
% Additional parameters - project specific
%-------------------------------------------------------------------------------%
numArrays                = 2; %number of arrays on each slide.
typeFlag                 = 'median'; % uses the median over replicates of an antigen. Can also be 'mean'
colorFlags               = {'635'};% <> can be changed to one or the other
colorTags                = {'IgG'}; %which dye was used for which Ab type.
switchIDs_and_Names_Flag = 0; % old flag for problem with Scrips GPR files, deprecated
commentFlag              = 1; % slide to sample mapping contain a comment column. Current default, used for backward compatbility

%handle to function that can process the GPR file. Each project will have a specific function.
processGPRfunctionHandle = @processIDRI_GPRstruct_Multiplex; 

%-------------------------------------------------------------------------------%

for i=1:length(experimentDates)%highlighted by zwm 07112014 i=4:4%length(experimentDates)
  mappingFilename = [loadPath,experimentDates{i},filesep,expPrefix{i},'SlideToSampleMappings.txt'];
  currLoadPath    = [loadPath,experimentDates{i},filesep,'GPR',filesep];
  filePrefix      = [expPrefix{i},'slide_']; % for legacy projects this may need to be modified

  % name of file in which converted data will be stored for future use
  arrayDataFilename = [savePath,filesep,experimentDates{i},filesep,expPrefix{i},'arrayData_',typeFlag,'.mat']; 
   
  [arrayData] = loadMultiplexArrayData(mappingFilename,filePrefix,fileSuffix,currLoadPath,numArrays,processGPRfunctionHandle,typeFlag,...
				       colorFlags,switchIDs_and_Names_Flag,commentFlag)
  save(arrayDataFilename,'arrayData');
end

