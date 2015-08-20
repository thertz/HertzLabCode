

% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

% To create new project files save as changing description, then search < > filling in with the appropriate info.
%
clear all; 
close all;

ZachFlag        = 0;
clusterPlotFlag = 1;
simplePlotFlag  = 1;
printFlag       = 0;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'KSHV';

if(ZachFlag)
 rootPath = [filesep,filesep,'fhdata.fhcrc.org',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,projectName,filesep];
  loadPath = savePath;
else % Tomer path's
  %rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  rootPath = [filesep,'Volumes',filesep,'Data',filesep,'HertzLab',filesep]; % main rootPath under which data is located  %local BGU path
  savePath = [rootPath,'ArrayData',filesep,projectName,filesep];
  loadPath = savePath;
end
%multiple dates can be entered, but must have a corresponding prefix
experimentDates = {'12_10_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'KSHV_12_10_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

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
expGroups                = {'pos','low_risk','NC'};% groups listed in sample to slide mapping
yLims                    = [0 25000];% max hight of y in graphs
minThreshold             = 10000; %minimal threshold for peak responses, used by findPeaks.
bgThreshold              = 1000; % threshold by which to flag antigens as ones that have high background and should
                                 % be removed from analysis (typically using a BSA negative control).
% Unsupervised clustering parameters:
numClusters          = 4; % number of clusters in the data (or expected number)
minResponseThreshold = 2000; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise
pathogenNames        = {'KSHV','Neg'}; %Neg18 is the E coli lysate alone.
numRepeats           = 3; %number ofreplicates of each antigen printed on the array.
figInd = 1;

% E. Coli lysate noramlization parameters:
EC_labels    = {'null', 'Dnull'};   % name of antigen label of the EC control spots
method      = 'max';  % summary statistic used to compute background, can be 'median' 'mean' or max. Note: Assumes that there are multiple NULL control spots on the chip
                      % at the same concentration, so if there are ones with different concentrations, will not work as expected. 
                      % Currently uses 'max' which bypasses the problem since we are using the highest concentration of the antigens (K_1_4)
EC_antigens = {'K','DK'}; % note that Null (the control antigen) is not listed here.

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

  
  % remove E.coli Lysate background noise by subtracting the response to E.Coli lysate for each individual from all
  % spots that were in this buffer. This is done in a sample specific manner, i.e we remove responses of bg 
  % using the control spot responses of that specific sample 
  [arrayData] = normalizeEC_lysateBGresponsesByBlock(arrayData, EC_labels, EC_antigens);
  
  slideNames            = [slideNames arrayData.slideNames];
  ptids                 = [ptids arrayData.ptids];
  groupNames            = [groupNames arrayData.groupNames];
  mPeptResponses        = [mPeptResponses ;[arrayData(1).responseMatrix_EC_subtracted{1}{1}]]; % NOTE USE OF SUBTRACTED ARRAY
  mPeptResponses(mPeptResponses < 0) = 0; %all negative responses zeroed out.

end

%-------------------------------------------------------------------------------%
% analysis 
%-------------------------------------------------------------------------------%


KSHVinds = [strmatch('K_',arrayData.antigenNames)];

figInd = 1;
if (clusterPlotFlag)   
  
  % cluster data using unserupvised correlation clustering and plot cluster dendrogram (tree) - need to define the
  % following parameters up-top: numClusters and minThreshold
  [distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses(:,KSHVinds),numClusters,minResponseThreshold);
  [corrMat] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,KSHVinds),minResponseThreshold);

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

  [figInd] = plotAllResponsesByGroup(mPeptResponses(:,KSHVinds),groupNames,expGroups,projectName,figInd,yLims,fontSize); %generic version

  % plot responses of eacn antigen comparing the specific groups
  %antigenNames = {'Orf 38','Orf 61', 'Orf 58','K5','K8.1','ORF 73','ORF 65','ORF 72','GFP','Null'}
  antigenNames = {'K_1', 'K_8', 'K_61','K_73'};
  
  %for i= 1:length(antigenNames)
  for i=KSHVinds'
    currInd = strmatch(arrayData.antigenNames{i},arrayData.antigenNames,'exact');
    figure(figInd);
    a = gca;
    set(a,'FontSize',16);  
    boxplot(mPeptResponses(:,currInd),groupNames);
    title(strrep(arrayData.antigenNames{i},'_',' '));
    figInd = figInd + 1;
  end
end

   %print figures to files:
if (printFlag)
    figPath = [savePath,'figs',filesep];
    for j=1:figInd
        figName = [figPath,projectName,'Figure_',num2str(j)];
        printFig(figure(j),figName,16,14);
    end
end
%}

