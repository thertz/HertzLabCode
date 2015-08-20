
% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

clear all;
close all;

ZachFlag        = 1;
simplePlotFlag  = 0;
clusterPlotFlag = 0;
%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'Dengue';

if(ZachFlag)
rootPath = [filesep,'\fhdata.fhcrc.org',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
savePath = [rootPath,'ArrayData',filesep,projectName,filesep];
loadPath = savePath;
else % Tomer path's
rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located
savePath = [rootPath,'ArrayData',filesep,projectName,filesep];
loadPath = savePath;
end

experimentDates = {'04_09_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'dengue_zm_04_09_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

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
fontSize = 16;
expGroups_h = {'NC_h','pre','post'};
expGroups_m = {'vac','NC_m','sublethal'};
yLims = [0 45000];
minThreshold             = 10000; %minimal threshold for peak responses, used by findPeaks.
% Unsupervised clustering parameters:
numClusters          = 3; % number of clusters in the data (or expected number)
minResponseThreshold = 1000; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise

%load Dengue peptide info
[pepData]   = readDengueNicaragua_PeptideData_z(rootPath);
Dengue1Inds = strmatch('DEN1',{pepData.strain});
Dengue2Inds = strmatch('DEN2',{pepData.strain});
Dengue3Inds = strmatch('DEN3',{pepData.strain});

Dengue1EnvInds = strmatch('DEN1|Env',{pepData.strain});
Dengue2EnvInds = strmatch('DEN2|Env',{pepData.strain});
Dengue3EnvInds = strmatch('DEN3|Env',{pepData.strain});
%-------------------------------------------------------------------------------%


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
  
 
  [peakResps] = findPeakResponses(mPeptResponses,minThreshold); %selects only responces over above set limit
  
 
  
  
end
%-------------------------------------------------------------------------------%
% Plotting-needs to be customized for each run
%-------------------------------------------------------------------------------%

figInd = 2;
if (simplePlotFlag)
    
 %plotting and save figures
  
  [figInd] = plotAllResponsesByGroup(mPeptResponses,groupNames,expGroups_h,projectName,figInd,yLims,fontSize);%for human samples
        printFig(2,[savePath,experimentDates{i},filesep,'FIG',filesep,'pre_post_h'],24,12);
  [figInd] = plotAllResponsesByGroup(mPeptResponses,groupNames,expGroups_m,projectName,figInd,yLims,fontSize);% for mouse
        printFig(3,[savePath,experimentDates{i},filesep,'FIG',filesep,'pre_post_m'],24,12);
  [figInd] = plotAllResponsesByPtid(mPeptResponses,ptids,'1264',projectName,figInd,yLims,fontSize);
        printFig(4,[savePath,experimentDates{i},filesep,'FIG',filesep,'1264'],24,12);
  [figInd] = plotAllResponsesByPtid(mPeptResponses,ptids,'1651',projectName,figInd,yLims,fontSize);
        printFig(5,[savePath,experimentDates{i},filesep,'FIG',filesep,'1651'],24,12);
  [figInd] = plotAllResponsesByPtid(mPeptResponses,ptids,'2008',projectName,figInd,yLims,fontSize);
        printFig(6,[savePath,experimentDates{i},filesep,'FIG',filesep,'2008'],24,12);

end  
       
       %--------------------------------------------------------------------%
if (clusterPlotFlag)    
  % cluster data using unserupvised correlation clustering and plot cluster dendrogram (tree) - need to define the
  % following parameters up-top: numClusters and minThreshold
  [distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses,numClusters,minResponseThreshold);
  [corrMat] = computeSampleCorrMatByResponseVectors(mPeptResponses,minResponseThreshold);
  
  % cluster with a sepcific set of antigens:
  [distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses(:,Dengue1Inds),numClusters,minResponseThreshold);
     
  % PLOT dendrogram:
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = ''; %dendrogram title
  [figInd] = plotResponseDendrogram(Zstruct,ptids,figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabels);
  [figInd] = plotSampleCorrmat(figInd,corrMat(I,I),ptids(I),[-1 1],fontSize);

end    
  
  
   %figure()
   %bar(vietInds(allContactInds),mPeptResponses([monoInds;baselineInds],vietInds(allContactInds))');
   %a = gca;
   %set(a,'FontSize',fontSize);
   %set(a,'yLim',[0 10000]);
   %legend(strrep(ptids([monoInds;baselineInds]),'_',' '));
  
%   figName = [savePath,experimentDates{i},filesep,expPrefix{i},monoNames{j},'_Responses'];
%   printFig(gcf,figName,20,16);

