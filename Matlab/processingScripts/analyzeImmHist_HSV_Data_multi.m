

% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

% To create new project files save as changing description, then search < > filling in with the appropriate info.
%
clear all;
close all;

ZachFlag        = 1;
clusterPlotFlag = 1;
simplePlotFlag  = 0;
printFlag       = 0;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'ImmuneHistory';

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
experimentDates = {'01_23_2014','02_05_2014','02_12_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'ImmuneHist_ZM_01_23_14_','ImmuneHist_ZM_02_05_14_','ImmuneHist_ZM_02_12_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

fileSuffix = ''; % deprcated option to have a suffix to the name - we no longer allow this in the current naming format.

%-------------------------------------------------------------------------------%
% Additional parameters - project specific
%-------------------------------------------------------------------------------%
numArrays                = 1; %number of arrays on each slide.
typeFlag                 = 'median'; % uses the median over replicates of an antigen. Can also be 'mean'
colorFlags               = {'635'};
colorTags                = {'IgG'};
switchIDs_and_Names_Flag = 0; % old flag for problem with Scrips GPR files, deprecated
commentFlag              = 1; % slide to sample mapping contain a comment column. Current default, used for backward compatbility
fontSize                 = 16;
expGroups                = {'HSV_1','HSV_2','HSV__1_2','HSV_neg','BSA'};% groups listed in sample to slide mapping
yLims                    = [0 35000];% max hight of y in graphs
minThreshold             = 10000; %minimal threshold for peak responses, used by findPeaks.
bgThreshold              = 1000; % threshold by which to flag antigens as ones that have high background and should
                                 % be removed from analysis (typically using a BSA negative control).
% Unsupervised clustering parameters:
numClusters          = 4; % number of clusters in the data (or expected number)
minResponseThreshold = 1000; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise
pathogenNames        = {'HSV_1_','HSV_2_','VZV','KSHV','Neg'}; %Neg18 is the E coli lysate alone.
numRepeats           = 3; %number ofreplicates of each antigen printed on the array.

%load Dengue peptide info- must be customized for new assay
%[pepData]   = <file must be created for each new type of slide> (rootPath);%not currently sure about this secion
%<Dengue1Inds = strmatch('DEN1',{pepData.strain});>
%<Dengue2Inds = strmatch('DEN2',{pepData.strain});>
%<Dengue3Inds = strmatch('DEN3',{pepData.strain});>

%<Dengue1EnvInds = strmatch('DEN1|Env',{pepData.strain});>
%<Dengue2EnvInds = strmatch('DEN2|Env',{pepData.strain});>
%<Dengue3EnvInds = strmatch('DEN3|Env',{pepData.strain});>
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

%should be fixed
         %HACK HACK HACK  - PROBLEM WITH NAMES OF ANTIGENS SEEM SLIKE BUT IN HSV_2 WHERE THERE IS NO HSV_2_11.
        % inds = strmatch('HSV_2_10',arrayData.seqGPRstructs{i}.Names);
        % arrayData.seqGPRstructs{1}.Names(inds(4:6)) = {'HSV_2_11'};
        % arrayData.antigenNames   = [arrayData.antigenNames(1:103); 'HSV_2_11';arrayData.antigenNames(104:end)];
        % arrayData.responseMatrix{1} = [arrayData.responseMatrix{1}(:,1:103) arrayData.responseMatrix{1}(:,103) arrayData.responseMatrix{1}(:,104:end)];
         %HACK HACK HACK  - PROBLEM WITH NAMES OF ANTIGENS SEEM SLIKE BUT IN HSV_2 WHERE THERE IS NO HSV_2_11.
  
  
  slideNames            = [slideNames arrayData.slideNames];
  ptids                 = [ptids arrayData.ptids];
  groupNames            = [groupNames arrayData.groupNames];
  mPeptResponses        = [mPeptResponses ;[arrayData(1).responseMatrix{1}]];

end

%get list of pathogens and their location on the mPeptResponse vectors:
[pathogenStruct] = getImmHistPathogenStructFromGPRstructNew(arrayData.seqGPRstructs{1},numRepeats,pathogenNames);


%remove E.coli Lysate background noise by subtracting the response to E.Coli lysate for each individual from all
%spots that were in this buffer
ECpathogens = {'Neg','HSV_1','HSV_2','KSHV','VZV'};

ECinds = [];
for i=1:length(ECpathogens)
  currInd = strmatch(ECpathogens{i},{pathogenStruct.name});  
  
  if(strcmp(ECpathogens{i},'Neg')) % look for Neg_18
    EC_lysateInd = strmatch('Neg_18',{pathogenStruct(currInd).antigens.name});
  end
  ECinds = [ECinds pathogenStruct(currInd).begInd:pathogenStruct(currInd).endInd];
end

  
for j=1:length(ptids)
  mPeptResponses(j,ECinds) = mPeptResponses(j,ECinds)-mPeptResponses(j,EC_lysateInd);
end

% subtract out negative control responses from all probes:
bgInds    = strmatch('HSV_neg',groupNames);
nonBGinds = setdiff([1:length(ptids)],bgInds);

bgResp = max(mPeptResponses(bgInds,:));

mPeptResponses(nonBGinds,:) = mPeptResponses(nonBGinds,:) - repmat(bgResp,length(nonBGinds),1);


%-------------------------------------------------------------------------------%
% analysis 
%-------------------------------------------------------------------------------%
%find peak responses to visually verify these on images:  
[peakResps] = findPeakResponses(mPeptResponses,minThreshold);
  

figInd = 1;
HSVinds = strmatch('HSV',arrayData.antigenNames);
KSHVinds = [strmatch('KSHV',arrayData.antigenNames)];
VZVinds = strmatch('VZV',arrayData.antigenNames);

if (clusterPlotFlag)   
  
  % cluster data using unserupvised correlation clustering and plot cluster dendrogram (tree) - need to define the
  % following parameters up-top: numClusters and minThreshold
  [distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses(:,HSVinds),numClusters,minResponseThreshold);
  [corrMat] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,HSVinds),minResponseThreshold);

  % cluster with a sepcific set of antigens:
  %[distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses(:,Dengue1Inds),numClusters,minResponseThreshold);
  
  % PLOT dendrogram:
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = ''; %dendrogram title
  [figInd] = plotResponseDendrogram(Zstruct,groupNames,figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabels);
  [figInd] = plotSampleCorrmat(figInd,corrMat(I,I),groupNames(I),[-1 1],fontSize);

end  
  %-------------------------------------------------------------------------------%  
  % plotting functions
  %-------------------------------------------------------------------------------%
if (simplePlotFlag)
  [figInd] = plotAllResponsesByGroup(mPeptResponses(:,HSVinds),groupNames,expGroups,'HSV-1 and HSV-2',figInd,yLims,fontSize);
  [figInd] = plotAllResponsesByGroup(mPeptResponses(:,VZVinds),groupNames,expGroups,'VZV',figInd,yLims,fontSize);
  %[figInd] = plotAllResponsesByPtid(mPeptResponses(:,HSVinds),ptids,'<>',projectName,figInd,yLims)

end


   %print figures to files:
if (printFlag)
    figPath = [savePath,filesep,'FIG',filesep];
    for j=1:figInd
        figName = [figPath,projectName,'Figure_',num2str(j)];
        printFig(figure(j),figName,30,24);
    end
end


