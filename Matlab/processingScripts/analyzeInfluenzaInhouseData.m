
% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

clear all;
close all;

ZachFlag = 0;


%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'Adjuvants';

if(ZachFlag)
 rootPath = [filesep,'\fhdata.fhcrc.org',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
else % Tomer path's
  rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
end
experimentDates = {'08_08_2013'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'Flu_inhouse_ZW_08_08_13_'}; % All slides and slideToSampleMapping file for each date have the exact same filePrefix

fileSuffix = ''; % deprcated option to have a suffix to the name - we no longer allow this in the current naming format.

%-------------------------------------------------------------------------------%
% Additional parameters - project specific
%-------------------------------------------------------------------------------%
numArrays                = 2; %number of arrays on each slide.
typeFlag                 = 'median'; % uses the median over replicates of an antigen. Can also be 'mean'
colorFlags               = {'635'};
colorTags                = {'IgG'}; %which dye was used for which Ab type.
switchIDs_and_Names_Flag = 0; % old flag for problem with Scrips GPR files, deprecated
commentFlag              = 1; % slide to sample mapping contain a comment column. Current default, used for backward compatbility
fontSize = 16;
expGroups = {'balbC','blk6','NC'};
yLims = [0 35000];




%handle to function that can process the GPR file. Each project will have a specific function.
processGPRfunctionHandle = @processFLU_HA_RapaGPRstruct_Multiplex; 

%-------------------------------------------------------------------------------%



%error here, need to change readHA_RapaPeptideData




% HA peptide info
[origArrayPepData] = readHA_RapaPeptideData(rootPath);
%peptides in idices 1:96 in vectors have the following tube labels (important for matching!)
vietInds = strmatch('Influenza A/VietNam/1203/2004|H5N1|HA ',{origArrayPepData.strain});
x31Inds  = strmatch('Influenza A/aichi/2/68 (X31 recomb)|H3N2|HA',{origArrayPepData.strain});
PR8Inds  = strmatch('Influenza A PR8/ThomasLab|HA|H1N1',{origArrayPepData.strain});

figInd = 1;
for i=1:length(experimentDates)
  mappingFilename = [loadPath,experimentDates{i},filesep,expPrefix{i},'SlideToSampleMappings.txt'];
  currLoadPath    = [loadPath,experimentDates{i},filesep,'GPR',filesep];
  filePrefix      = [expPrefix{i},'slide_'];

  % name of file in which converted data will be stored for future use
  arrayDataFilename = [savePath,experimentDates{i},filesep,expPrefix{i},'arrayData_',typeFlag,'.mat']; 

  load(arrayDataFilename);
 
  slideNames = [arrayData.slideNames];
  ptids      = [arrayData.ptids];
  groupNames = [arrayData.groupNames];
  mPeptResponses      = [arrayData(1).responseMatrix{1}];

  % plot all responses by group - limited only to the viet1203 inds
  [figInd] = plotAllResponsesByGroup(mPeptResponses(:,vietInds),groupNames,expGroups,'Adjuvants',figInd,yLims,fontSize); 
  
  % plot all responses by group after background subtraction
  NCinds = strmatch('NC',groupNames);
  mPeptResponsesBgSubtracted = mPeptResponses -repmat(max(mPeptResponses(NCinds,:)),24,1);

  % put back the NC responses unsubtracted!
  mPeptResponsesBgSubtracted(NCinds,:) = mPeptResponses(NCinds,:);
  [figInd] = plotAllResponsesByGroup(mPeptResponsesBgSubtracted(:,vietInds),groupNames,expGroups,'Adjuvants',figInd,yLims,fontSize);

  [figInd] = plotAllResponsesByPtid(mPeptResponses(:,vietInds),ptids,{'balbc_10','balbc_5'},'Adjuvants',figInd,yLims)

  
  %print figures to files:
  figPath = [savePath,experimentDates{i},filesep,'figures',filesep];
  for j=1:3
    figName = [figPath,projectName,'Figure_',num2str(j)];
    printFig(figure(j),figName);
  end
end

