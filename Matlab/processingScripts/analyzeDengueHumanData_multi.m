
% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.
%
% for Dengue humans samples from the Harris Lab - updated to work on new data generated on in-house printed slides Dec. 2014


clear all;
close all;

ZachFlag        = 0;
plotFlag        = 1;
printFlag       = 0;
%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'Dengue';

if(ZachFlag)
  rootPath = [filesep,'\fhdata.fhcrc.org',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
else % Tomer path's
  %rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located
  %rootPath = [filesep,'Volumes',filesep,'Data',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  rootPath = ['~/Dropbox',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
end
savePath = [rootPath,'ArrayData',filesep,projectName,filesep];
loadPath = savePath;

% 04_07_2014 - experiment comparing serum to purified IgG using Raining tips.
%experimentDates = {'04_07_2014','06_27_2014','08_28_2014','04_03_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
%expPrefix    = {'dengue_zm_04_07_14_','dengue_zm_06_27_14_','dengue_zm_08_28_14_','dengue_zm_04_03_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

% new samples run on new in-house printed dengue slides;
experimentDates = {'12_01_2014','12_03_2014','12_05_2014','12_09_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'dengue_zm_12_01_14_','dengue_zm_12_03_14_','dengue_zm_12_05_14_','dengue_zm_12_09_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix
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
if(printFlag)
  fontSize = 8;
else
  fontSize = 16;
end

%expGroups = {'Human_Convalescent','Human_Acute','Human_post','Human_pre', 'Human_12Month','Human_Naive'}; % for new runs on in-house slides Dec. 2014
expGroups = {'Human_Acute', 'Human_Convalescent','Human_12Month','Human_Naive'}; % for new runs on in-house slides Dec. 2014

%yLims = [log2(1024) log2(65536)];
%minThreshold             = log2(1024); %minimal threshold for peak responses, used by findPeaks.
yLims        = [0 65000];
minThreshold = 1000;


% Unsupervised clustering parameters:
numClusters          = 6; % number of clusters in the data (or expected number)
%minResponseThreshold = log2(2048); %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise
minResponseThreshold = 2000;

%-------------------------------------------------------------------------------%
% Stat comparison parameters: (also project specific)
%-------------------------------------------------------------------------------%
% Filters used for treatment-blinded filtering of responses to single antigens used for statistical analysis.
% currently only one filter is implemented which is % responders within the entire treatment set (ignoring all controls) 
filterNames = {'PercentResponders'};
filterThresholds = {0.2};
treatmentGroups = {'Human_Acute','Human_Convalescent'}; % Names of the two groups that are to be compared statistically below (currently supports only pairwise comparisons)


%load Dengue peptide info
[pepData]   = readDengueNicaragua_PeptideData(rootPath);
Dengue1Inds = strmatch('DEN1',{pepData.strain});
Dengue2Inds = strmatch('DEN2',{pepData.strain});
Dengue3Inds = strmatch('DEN3',{pepData.strain});

Dengue1EnvInds = strmatch('DEN1|Env',{pepData.strain});
Dengue2EnvInds = strmatch('DEN2|Env',{pepData.strain});
Dengue3EnvInds = strmatch('DEN3|Env',{pepData.strain});
%-------------------------------------------------------------------------------%
%defining elements outside of loop to add loop data too

slideNames      = {};
ptids           = {};
groupNames      = {};
mPeptResponses  = [];
expDates        = []
%-------------------------------------------------------------------------------%
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
  expDates              = [expDates ones(1,length(arrayData.ptids))*i];
 
end

% first clean data to include only human samples, as some runs included both human/mouse
humanInds = strmatch('Human',groupNames);

slideNames     = slideNames(humanInds);
ptids          = ptids(humanInds);
groupNames     = groupNames(humanInds);
mPeptResponses = mPeptResponses(humanInds,:);

% now recompute indices of neg controls after removing all mouse/mAb data:
NegControlInds = strmatch('Human_Naive',groupNames); % 
infectedInds   = setdiff([1:length(ptids)],NegControlInds);

% small HACK HACK HACK to correct one bad gropuName label:
badInd = strmatch('Human_Human_Acute',groupNames);
groupNames{badInd} = 'Human_Acute';

% normalize background responses - currently using ajduvant with OVA:
%NegControlInds = strmatch('Naive',ptids); % Older version for previous dates...

bgResponses    = median(mPeptResponses(NegControlInds,:));
mPeptResponsesNormalized = mPeptResponses; 
mPeptResponsesNormalized(infectedInds,:) = mPeptResponses(infectedInds,:)-repmat(bgResponses,length(infectedInds),1);

%another optinal way to compute background when two channels are available;
%[naiveResponses] = computeBgResposnesFromArrayDataStruct(arrayData,'OVA_MF59');
  
% clean data and take log10 transform:
%mPeptResponses(mPeptResponses<1024) = 1024;
%mPeptResponses = log2(mPeptResponses);

%mPeptResponsesNormalized(mPeptResponsesNormalized<1024) = 1024;
%mPeptResponsesNormalized = log2(mPeptResponsesNormalized);

mPeptResponses(mPeptResponses < 0)                     = 0;
mPeptResponsesNormalized(mPeptResponsesNormalized < 0) = 0;


%-------------------------------------------------------------------------------%
% Compute pairwise statistics comparing treatment groups of interest
% 
% Note: We are comparing normalized (background subtracted) responses.
%-------------------------------------------------------------------------------%
inds1 = strmatch(treatmentGroups{1},groupNames)';
inds2 = strmatch(treatmentGroups{2},groupNames)';
if(isempty(inds1) || isempty(inds2))
  disp(sprintf('Treatment groups for stat comparison not found'))
  statFlag = 0;
else
  statFlag = 1;
  respMat = mPeptResponsesNormalized([inds1,inds2],:);
  treatmentLabels  = [zeros(1,length(inds1))  ones(1,length(inds2))];

  % compute label blinded filtering of antigens across the array - currenlty using %responders - set at 20%
  [antigenFilter] = computeRespMatLabelBlindedFilters(respMat,filterNames,filterThresholds,minResponseThreshold);
  filteredAntigenInds = find(antigenFilter); %list of indices to antigens in which more than % repsonderes responded (treatment blinded)

  % compute analysis using responses to all 3 Dengue strains:
  [groupStatsOverall] = compareTreatmentGroups(respMat,treatmentLabels,minResponseThreshold);

  %limit analysis to Dengue2 (vaccine strain):
  [groupStatsDengue2] = compareTreatmentGroups(respMat(:,Dengue2Inds),treatmentLabels,minResponseThreshold);


  % OUTPUT results to command line:
  disp(sprintf('*----------------------------------------------------------------------------*'))
  disp(sprintf('Comparing overall responses across the entire Dengue array'))
  disp(sprintf('Ranksum p-value: %0.5g\n',groupStatsOverall.totResponses_pValue))
  disp(sprintf('*----------------------------------------------------------------------------*'))

  disp(sprintf('*----------------------------------------------------------------------------*'))
  disp(sprintf('Comparing Dengue2 responses'))
  disp(sprintf('Ranksum p-value: %0.5g\n',groupStatsDengue2.totResponses_pValue))
  disp(sprintf('*----------------------------------------------------------------------------*'))

  disp(sprintf('*----------------------------------------------------------------------------*'))
  disp(sprintf('Ranksum p-values for filtered antigens:'))
  for i=1:length(filteredAntigenInds)
    currInd = filteredAntigenInds(i);
    disp(sprintf('Peptide #:%d, Strain|protein: %s, startPosition: %d, Fisher p-value: %0.5g',...
      currInd,pepData(currInd).strain, pepData(currInd).begInd, groupStatsOverall.singleAntigen_pValuesFisher(currInd)));
  end
  disp(sprintf('*----------------------------------------------------------------------------*'))
end


% unsupervised clustering:
[distMatDen2,ZstructDen2,clusterLabelsDen2] = clusterSamplesByResponseVectors(mPeptResponses(:,Dengue2Inds),numClusters,minResponseThreshold);
[corrMatDen2] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,Dengue2Inds),minResponseThreshold);      

%-------------------------------------------------------------------------------%
% Plotting-needs to be customized for each run
%-------------------------------------------------------------------------------%
figInd = 1;
if (plotFlag)
    
 %plotting and save figures
  [figInd] = plotAllResponsesByGroup(mPeptResponsesNormalized,groupNames,expGroups,projectName,figInd,yLims,fontSize);

  % PLOT dendrogram: - currently ploting only overall clustering results
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = ''; %dendrogram title
  [figInd] = plotResponseDendrogram(ZstructDen2,groupNames,figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabelsDen2);
  [figInd] = plotSampleCorrmat(figInd,corrMatDen2(I,I),groupNames(I),[0.5 1],fontSize);

  % plot responses by clusters
  [figInd] = plotAllResponsesByClusters(mPeptResponses(:,Dengue2Inds),clusterLabelsDen2,'Dengue 2',figInd,yLims,fontSize)

  
  % plot boxplots of overall responses comparing the two treatment groups
  [figInd] = plotGroupResponsesBoxPlots(groupStatsOverall.totResponses,treatmentLabels,treatmentGroups,...
    groupStatsOverall.totResponses_pValue,figInd,fontSize);

  % plot individual significant antigens:
  if(statFlag)
    [figInd,sigInds] = plotIndividualResponsesByGroup(respMat,treatmentLabels,treatmentGroups,...
      groupStatsOverall.singleAntigen_pValuesFisher,0.05,Dengue1Inds,'Dengue1',pepData,figInd,fontSize);

    [figInd,sigInds] = plotIndividualResponsesByGroup(respMat,treatmentLabels,treatmentGroups,...
      groupStatsOverall.singleAntigen_pValuesFisher,0.05,Dengue2Inds,'Dengue2',pepData,figInd,fontSize);

    [figInd,sigInds] = plotIndividualResponsesByGroup(respMat,treatmentLabels,treatmentGroups,...
      groupStatsOverall.singleAntigen_pValuesFisher,0.05,Dengue3Inds,'Dengue3',pepData,figInd,fontSize);
  end
  % plot responses of antigens of interest (not necessarily significant)
  %antigenInd = 330;
  %[figInd] = plotIndividualResponsesForSingleAntigen(respMat(:,antigenInd),antigenInd,treatmentLabels,...
  %  treatmentGroups,pepData(antigenInd),figInd,fontSize);

  % plot responses to the 4 peptides identified by the mouse data and ordere from CPC (LAST PEPTIDES ON THIS ARRAY NOT ANNOTATED IN PEPDATA YET)
  mousePepInds = [345,347,330, 356, 588:591]; % original peptides and last 4 peptides on the array
  antigenNames = {'D345','D347', 'D330 ', 'D356','D345 - New','D347 - New', 'D330 - New', 'D356 - New' };
  [figInd] = plotIndividualBarPeptideResponsesByGroup(mPeptResponses(:,mousePepInds),groupNames,expGroups,antigenNames,yLims,figInd,fontSize)

  % plot boxplot comparisons of response in each group to each peptide:
  allGroupInds    = [];
  allGroupLabels  = [];
  for i=1:length(expGroups)
    currGroupInds = strmatch(expGroups{i},groupNames);
    allGroupInds = [allGroupInds ; currGroupInds];
    allGroupLabels = [allGroupLabels ones(1,length(currGroupInds))*(i-1)];
  end

  for i=1:length(mousePepInds)
    [figInd] = plotGroupResponsesBoxPlots(mPeptResponses(allGroupInds,mousePepInds(i)),allGroupLabels,strrep(expGroups,'Human_',''),...
      groupStatsOverall.singleAntigen_pValuesFisher(mousePepInds(i)),figInd,fontSize,antigenNames{i});
      set(gcf,'Position',[300,400,800,700]);
  end

end 
 

if (printFlag)
  figPath = [savePath,'figs',filesep,'Human',filesep];
  for j=1:figInd
    figName = [figPath,projectName,'Figure_',num2str(j)];
    [figHandle] = printFigure(j,figName,{'eps','png'});
    close(figure(j));
  end
end
