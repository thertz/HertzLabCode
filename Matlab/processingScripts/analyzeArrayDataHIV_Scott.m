

% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

% To create new project files save as changing description, then search < > filling in with the appropriate info.
%
clear all;
close all;

ZachFlag  = 1;
clusterPlotFlag = 0;
simplePlotFlag  = 0;
printFlag = 0;
zBar    = 0;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'HIV_Scott';

if(ZachFlag)
 rootPath = [filesep,filesep,'fhdata.fhcrc.org',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,projectName,filesep];
  loadPath = savePath;
else % Tomer path's
  rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,projectName,filesep];
  loadPath = savePath;
end
%multiple dates can be entered, but must have a corresponding prefix
experimentDates = {'08_14_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'HIV_08_14_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

fileSuffix = ''; % deprcated option to have a suffix to the name - we no longer allow this in the current naming format.

%-------------------------------------------------------------------------------%
% Additional parameters - project specific
%-------------------------------------------------------------------------------%
numArrays                = 1; %number of arrays on each slide.
typeFlag                 = 'mean'; % uses the median over replicates of an antigen. Can also be 'mean'
colorFlags               = {'635'};
colorTags                = {'IgG'};
switchIDs_and_Names_Flag = 0; % old flag for problem with Scrips GPR files, deprecated
commentFlag              = 1; % slide to sample mapping contain a comment column. Current default, used for backward compatbility
fontSize                 = 16;
expGroups                = {};% groups listed in sample to slide mapping
yLims                    = [0 65000];% max hight of y in graphs
minThreshold             = 10000; %minimal threshold for peak responses, used by findPeaks.
bgThreshold              = 5000; % threshold by which to flag antigens as ones that have high background and should
rsGroups                 ={ 'NC','Elite_controlers', 'Slow_progressor','Rapid_progressor'}  ; 
seroNormGroups           ={'NC',  'Normal_progressor','Seroconverter'};
multiGroups              = {'Seroconverter','Slow_progressor','Normal_progressor','Rapid_progressor'};

% Unsupervised clustering parameters:
numClusters          = 3; % number of clusters in the data (or expected number)
minResponseThreshold = 1000; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise

%load Dengue peptide info- must be customized for new assay
%{
[pepData]   = <file must be created for each new type of slide> (rootPath);%not currently sure about this secion
<Dengue1Inds = strmatch('DEN1',{pepData.strain});>
<Dengue2Inds = strmatch('DEN2',{pepData.strain});>
<Dengue3Inds = strmatch('DEN3',{pepData.strain});>

<Dengue1EnvInds = strmatch('DEN1|Env',{pepData.strain});>
<Dengue2EnvInds = strmatch('DEN2|Env',{pepData.strain});>
<Dengue3EnvInds = strmatch('DEN3|Env',{pepData.strain});>
%}
%-------------------------------------------------------------------------------%
%defining elements outside of loop to add loop data too
slideNames      ={};
ptids           ={};
groupNames      ={};
mPeptResponses  =[];
%-------------------------------------------------------------------------------


for i=1:length(experimentDates)
  mappingFilename = [loadPath,experimentDates{i},filesep,expPrefix{i},'SlideToSampleMappings.txt'];
  currLoadPath    = [loadPath,experimentDates{i},filesep,'GPR',filesep];
  filePrefix      = [expPrefix{i},'slide_'];

  % name of file in which converted data will be stored for future use
  arrayDataFilename = [savePath,experimentDates{i},filesep,expPrefix{i},'arrayData_',typeFlag,'.mat']; 

  load(arrayDataFilename);
 
  slideNames            = [slideNames arrayData.slideNames];
  ptids                 = [ptids arrayData.ptids];
  groupNames            = [groupNames arrayData.groupNames];
  mPeptResponses        = [mPeptResponses ;[arrayData(1).responseMatrix{1}]];

end
figInd = 1;
  %-------------------------------------------------------------------------------%
  % analysis 
  %-------------------------------------------------------------------------------%
  %find peak responses to visually verify these on images:
  
  [peakResps] = findPeakResponses(mPeptResponses,minThreshold);
if (clusterPlotFlag)   
  
  % cluster data using unserupvised correlation clustering and plot cluster dendrogram (tree) - need to define the
  % following parameters up-top: numClusters and minThreshold
  [distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses,numClusters,minResponseThreshold);
  [corrMat] = computeSampleCorrMatByResponseVectors(mPeptResponses,minResponseThreshold);

  % cluster with a sepcific set of antigens:
  %[distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses(:,Dengue1Inds),numClusters,minResponseThreshold);
  
  % PLOT dendrogram:
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = ''; %dendrogram title
  [figInd] = plotResponseDendrogram(Zstruct,ptids,figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabels);
  [figInd] = plotSampleCorrmat(figInd,corrMat(I,I),ptids(I),[-1 1],fontSize);

end  
  %-------------------------------------------------------------------------------%  
  % plotting functions
  %-------------------------------------------------------------------------------%
if (simplePlotFlag)
  [figInd] = plotAllResponsesByGroup(mPeptResponses,groupNames,expGroups,projectName,figInd,yLims,fontSize);%
  [figInd] = plotAllResponsesByPtid(mPeptResponses,ptids,'<>',projectName,figInd,yLims)

end
if (zBar)
    [figInd] = barPlotGroups(mPeptResponses,groupNames,rsGroups,'HIV Antigens',figInd,yLims, fontSize)
    [figInd] = barPlotGroups(mPeptResponses,groupNames,seroNormGroups ,'HIV Antigens',figInd,yLims, fontSize)
    [figInd] = barPlotGroups(mPeptResponses,groupNames,multiGroups ,'HIV Antigens',figInd,yLims, fontSize)
end
   %print figures to files:
if (printFlag)
    figPath = [savePath,filesep,'FIG',filesep];
    for j=1:figInd
        figName = [figPath,projectName,'Figure_',num2str(j)];
        printFig(figure(j),figName,30,24);
    end
end


