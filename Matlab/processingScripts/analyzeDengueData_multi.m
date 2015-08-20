
% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

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
  rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located
end
savePath = [rootPath,'ArrayData',filesep,projectName,filesep];
loadPath = savePath;

% mAb 1H7.4 experiments with different dilutions:
experimentDates = {'07_15_2014','07_16_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'dengue_zm_07_15_14_','dengue_zm_07_16_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

% mouse experiments, includin older (first two dates) and newer samples
% dates 1-2 are from first run on older slides. date 3 is the new run on new slides with new samples. (they cannot be analyzed together)
%experimentDates = {'08_21_2013','12_23_2013','04_09_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
%expPrefix    = {'dengue_zm_08_21_13_','dengue_zm_12_23_13_','dengue_zm_04_09_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

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
expGroupsMouse = {'MF59','Plo46','Alum','OVA_MF59','NMS','BSA'};

yLims = [log2(1024) log2(65536)];
minThreshold             = log2(1024); %minimal threshold for peak responses, used by findPeaks.
% Unsupervised clustering parameters:
numClusters          = 4; % number of clusters in the data (or expected number)
minResponseThreshold = log2(4096); %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise

%-------------------------------------------------------------------------------%
% Stat comparison parameters: (also project specific)
%-------------------------------------------------------------------------------%
% Filters used for treatment-blinded filtering of responses to single antigens used for statistical analysis.
% currently only one filter is implemented which is % responders within the entire treatment set (ignoring all controls) 
filterNames = {'PercentResponders'};
filterThresholds = {0.2};
treatmentGroups = {'Alum','MF59'}; % Names of the two groups that are to be compared statistically below (currently supports only pairwise comparisons)

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
slideNames     = {};
ptids          = {};
groupNames     = {};
mPeptResponses = [];
arrayData      =[];
%-------------------------------------------------------------------------------%
for i=1:length(experimentDates)
  mappingFilename = [loadPath,experimentDates{i},filesep,expPrefix{i},'SlideToSampleMappings.txt'];
  currLoadPath    = [loadPath,experimentDates{i},filesep,'GPR',filesep];
  filePrefix      = [expPrefix{i},'slide_'];

  % name of file in which converted data will be stored for future use
  arrayDataFilename = [savePath,experimentDates{i},filesep,expPrefix{i},'arrayData_',typeFlag,'.mat']; 

  currArrayData = load(arrayDataFilename);

  %small hack to change the Adjuvant only (negative control) ptids to another name:
  inds = strmatch('MF59_plus_OVA',currArrayData.arrayData.ptids);
  currArrayData.arrayData.ptids(inds) = {'OVA_MF59'};

  slideNames            = [slideNames currArrayData.arrayData.slideNames];
  ptids                 = [ptids currArrayData.arrayData.ptids];
  groupNames            = [groupNames currArrayData.arrayData.groupNames];
  mPeptResponses        = [mPeptResponses ;[currArrayData.arrayData(1).responseMatrix{1}]];
 
  arrayData = [arrayData currArrayData.arrayData];
  %[peakResps] = findPeakResponses(mPeptResponses,minThreshold); %selects only responces over above set limit
  
end

% normalize background responses - currently using ajduvant with OVA:
NegControlInds = strmatch('OVA_MF59',ptids); % Adjuvant and Ova negative control
bgResponses    = mean(mPeptResponses(NegControlInds,:));
mPeptResponsesNormalized = mPeptResponses-repmat(bgResponses,size(mPeptResponses,1),1);

%another optinal way to compute background when two channels are available;
%[naiveResponses] = computeBgResposnesFromArrayDataStruct(arrayData,'OVA_MF59');
  
% clean data and take log10 transform:
mPeptResponses(mPeptResponses<1024) = 1024;
mPeptResponses = log2(mPeptResponses);

mPeptResponsesNormalized(mPeptResponsesNormalized<1024) = 1024;
mPeptResponsesNormalized = log2(mPeptResponsesNormalized);

%-------------------------------------------------------------------------------%
% Compute pairwise statistics comparing treatment groups of interest
% 
% Note: We are comparing normalized (background subtracted) responses.
%-------------------------------------------------------------------------------%
inds1 = strmatch(treatmentGroups{1},groupNames)';
inds2 = strmatch(treatmentGroups{2},groupNames)';
if(isempty(inds1) || isempty(inds2))
  disp(sprintf('Treatment groups for stat comparison not found'))
  StatFlag = 0;
else
  StatFlag = 1;
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

%----------------------------------------------------------------------------------------------------------------%
% cluster data using unserupvised correlation clustering and plot cluster dendrogram (tree) - need to define the
% following parameters up-top: numClusters and minThreshold
%----------------------------------------------------------------------------------------------------------------%
% using all antigens
[distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses,numClusters,minResponseThreshold);
[corrMat] = computeSampleCorrMatByResponseVectors(mPeptResponses,minResponseThreshold);
  
% cluster with a sepcific set of antigens: (Dengue 2)
[distMat1,Zstruct1,clusterLabels1] = clusterSamplesByResponseVectors(mPeptResponses(:,Dengue2Inds),numClusters,minResponseThreshold);
[corrMat1] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,Dengue2Inds),minResponseThreshold);


%-------------------------------------------------------------------------------%
% Plotting-needs to be customized for each run
%-------------------------------------------------------------------------------%
figInd = 1;
if (plotFlag)
    
 %plotting and save figures
  [figInd] = plotAllResponsesByGroup(mPeptResponses,ptids,expGroupsMouse,projectName,figInd,yLims,fontSize);% for mouse
  [figInd] = plotAllResponsesByGroup(mPeptResponsesNormalized,ptids,{'MF59','Alum'},projectName,figInd,yLims,fontSize);% for mouse

      
  % PLOT dendrogram: - currently ploting only overall clustering results
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = ''; %dendrogram title
  [figInd] = plotResponseDendrogram(Zstruct,ptids,figInd,fontSize,bwFlag,titleStr);
  %[figInd] = plotResponseDendrogram(Zstruct1,ptids,figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabels);
  [figInd] = plotSampleCorrmat(figInd,corrMat(I,I),ptids(I),[0 1],fontSize);
  %[Y I] = sort(clusterLabels1);
  %[figInd] = plotSampleCorrmat(figInd,corrMat1(I,I),ptids(I),[0 1],fontSize);

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

end 
 
%print figures to files:
if (printFlag)
  figPath = [savePath,'figs',filesep];
  for j=1:figInd
    figName = [figPath,projectName,'Figure_',num2str(j)];
    printFig(figure(j),figName,30,24);
  end
end
