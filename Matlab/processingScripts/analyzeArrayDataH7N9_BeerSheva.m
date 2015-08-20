

% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

% To create new project files save as changing description, then search < > filling in with the appropriate info.
%
clear all;
close all;

ZachFlag        = 1;
clusterPlotFlag = 0;
plotFlag        = 0;
printFlag       = 0;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'H7N9';

if(ZachFlag)
  rootPath = ['C:\Users\friedmal\Documents\Tomer Hertz\experiments\'];
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
else % Tomer path's
  %rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  rootPath = [filesep,'Volumes',filesep,'Data',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
end
%multiple dates can be entered, but must have a corresponding prefix
experimentDates = {'06_03_2015'};%{'08_21_2014','08_22_2014','08_25_2014','09_04_2014','09_05_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'H7N9_06_03_15_'};%{'H7N9_08_21_14_','H7N9_08_22_14_','H7N9_08_25_14_','H7N9_09_04_14_','H7N9_09_05_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix


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
fontSize                 = 12;

% Q: - What are OB_Vac_post and WT_Vac_post etc?
adjuvantLabels    = {'Vac','AS03','MF59','PBS'}; % types of adjuvants used on both Obest and WT mice
expGroupPrefixes = {'WT_pre_','Ob_pre_','WT_post_','Ob_post_'};
    
yLims                    = [0 65000];% max height of y in graphs
yLimsSummary             = [0 32500]; %max height of y in summary stat graphs
minThreshold             = 10000; %minimal threshold for peak responses, used by findPeaks.
bgThreshold              = 5000; % threshold by which to flag antigens as ones that have high background and should
                                 % be removed from analysis (typically using a BSA negative control).
summaryMethod = 'median';  % summary stat used for group comparisons and for plotting gropu responses - and cluster responses

% Unsupervised clustering parameters:
numClusters          = 5; % number of clusters in the data (or expected number)
minResponseThreshold = 2000; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise

%-------------------------------------------------------------------------------%
% Stat comparison parameters: (also project specific)
%-------------------------------------------------------------------------------%
% Filters used for treatment-blinded filtering of responses to single antigens used for statistical analysis.
% currently only one filter is implemented which is % responders within the entire treatment set (ignoring all controls) 
filterNames = {'PercentResponders'};
filterThresholds = {0.2};
treatmentGroups = {'Ob_post_Vac','WT_post_Vac'}; % Names of the two groups that are to be compared statistically below (currently supports only pairwise comparisons)

% data Labels
maskLabels = {'Cal_HA_Inds','Cal_NA_Inds','Sha_HA_Inds','Sha_NA_Inds'};
treatmentLabels = {'WT-pre-vac', 'obese-pre-vac','wt-post-vac','obese-post-vac'};

%-------------------------------------------------------------------------------%
%defining elements outside of loop to add loop data too
slideNames     = {};
ptids          = {};
groupNames     = {};
mPeptResponses = [];
expDates       = [];
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
  expDates              = [expDates ones(1,length(arrayData.ptids))*i];
end
figInd = 1;
%{
%-----------------------------------------------------------------%
% Background subtraction
%-----------------------------------------------------------------%
%remove negative control responses as measured by BSA alone:
BSAinds     = strmatch('BSA',ptids)';
nonBSAinds  = setdiff([1:length(ptids)],BSAinds);
if(length(BSAinds)>1)
  bgResponses = max(mPeptResponses(BSAinds,:));
else
  bgResponses = mPeptResponses(BSAinds,:);
end
  
% split out BSA responses - to allow label consistency with the FIM12struct
BSAresponses   = mPeptResponses(BSAinds,:);

slideNames = slideNames(nonBSAinds);
ptids = ptids(nonBSAinds);
groupNames = groupNames(nonBSAinds);

mPeptResponses = mPeptResponses(nonBSAinds,:);
mPeptResponses = mPeptResponses - repmat(bgResponses,length(nonBSAinds),1);

% clean data and take log2 transform:
mPeptResponses(mPeptResponses<0) = 0;
%mPeptResponses = log2(mPeptResponses);


% Organize antigen mask sets
antigenNames = arrayData(1).antigenNames;
strainNames  = {'Cal','Sha'};    
%index sets for the different strains/proteins:
Cal_HA_Inds = strmatch('Cal_ha',antigenNames);
Cal_NA_Inds = strmatch('Cal_na',antigenNames);
Sha_HA_Inds = strmatch('SHA_ha',antigenNames);
Sha_NA_Inds = strmatch('SHA_na',antigenNames);


% Struct use to encode all antigen masks to allwo streamlining code for iterating over mask sets
antigenMasksStruct(1).name = maskLabels{1};
antigenMasksStruct(1).inds = Cal_HA_Inds;

antigenMasksStruct(2).name = maskLabels{2};
antigenMasksStruct(2).inds = Cal_NA_Inds;

antigenMasksStruct(3).name = maskLabels{3};
antigenMasksStruct(3).inds = Sha_HA_Inds;

antigenMasksStruct(4).name = maskLabels{4};
antigenMasksStruct(4).inds = Sha_NA_Inds;

CalInds = [Cal_HA_Inds;Cal_NA_Inds];
ShaInds = [Sha_HA_Inds; Sha_NA_Inds];

for i=1:length(adjuvantLabels)
  expGroupStruct(i).adjuvant   = adjuvantLabels{i};
  expGroupStruct(i).groupNames = {};
  expGroupStruct(i).inds       = [];
  for j=1:length(expGroupPrefixes)
    currGroup                      = [expGroupPrefixes{j},adjuvantLabels{i}];
    expGroupStruct(i).groupNames   = [expGroupStruct(i).groupNames currGroup];
    currInds                       = strmatch(lower(currGroup),lower(groupNames));
    expGroupStruct(i).inds         = [expGroupStruct(i).inds ;currInds];
    expGroupStruct(i).groupInds{j} = currInds; % save by group for different slicing.
    expGroupStruct(i).labels(currInds) = j; % label vec for boxplot comparisons
  end
end



  %-------------------------------------------------------------------------------%
  % analysis 
  %-------------------------------------------------------------------------------%

  % cluster just post boost responses
  postInds = [strmatch('Ob_post',groupNames) ; strmatch('WT_post',groupNames)];

  [distMatCal,ZstructCal,clusterLabelsCal] = clusterSamplesByResponseVectors(mPeptResponses(postInds,CalInds),numClusters,minResponseThreshold,'spearman','complete');
  [corrMatCal] = computeSampleCorrMatByResponseVectors(mPeptResponses(postInds,CalInds),minResponseThreshold);

  [distMatSha,ZstructSha,clusterLabelsSha] = clusterSamplesByResponseVectors(mPeptResponses(postInds,ShaInds),numClusters,minResponseThreshold,'spearman','complete');
  [corrMatSha] = computeSampleCorrMatByResponseVectors(mPeptResponses(postInds,ShaInds),minResponseThreshold);

  %now plot for each adjuvant separately, retaining only post boost samples again.
  for i=1:length(expGroupStruct)

    currInds = intersect(expGroupStruct(i).inds,postInds);

    [distMatGroups{i},ZstructGroups{i},clusterLabelsGroups{i}] = ...
       clusterSamplesByResponseVectors(mPeptResponses(currInds,ShaInds),numClusters,minResponseThreshold,'spearman','complete');
    [corrMatGroups{i}] = computeSampleCorrMatByResponseVectors(mPeptResponses(currInds,ShaInds),minResponseThreshold);
  end


if (clusterPlotFlag)  
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabelsSha);
  [figInd] = plotSampleCorrmat(figInd,corrMatSha(I,I),groupNames(postInds(I)),[0 1],fontSize);
  titleStr = 'Shanghai H7N9'; %dendrogram title
  title(titleStr);

  % PLOT dendrograms:
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = 'Shanghai H7N9'; %dendrogram title
  [figInd] = plotResponseDendrogram(ZstructSha,groupNames(postInds),figInd,fontSize-2,bwFlag,titleStr);

  for i=1:length(expGroupStruct)
    currInds = intersect(expGroupStruct(i).inds,postInds);
    [figInd] = plotResponseDendrogram(ZstructGroups{i},groupNames(currInds),figInd,fontSize-2,bwFlag,titleStr);
  end 

  % plot the summary (mean/median) response of each of the clusters induced by Shanghai :
  [figInd] = plotSummaryResponsesByClusters(mPeptResponses,clusterLabelsSha,titleStr,summaryMethod,figInd,yLims,fontSize);
  [figInd] = plotAllResponsesByClusters(mPeptResponses,clusterLabelsSha,titleStr,figInd,yLims,fontSize);
end


  %-------------------------------------------------------------------------------%  
  % plotting functions
  %-------------------------------------------------------------------------------%
if (plotFlag)

  % compare an adjuvant across different groups:
  for i=1:length(expGroupStruct)
    %[figInd] = plotAllResponsesByGroup(mPeptResponses,groupNames,[expGroupStruct(i).groupNames, 'NC'],projectName,figInd,yLims,fontSize);
    [figInd] = plotSummaryResponsesByGroup(mPeptResponses,groupNames,[expGroupStruct(i).groupNames, 'NC'],projectName,summaryMethod,figInd,yLimsSummary,fontSize);
  end

  % plot responses comparing different adjuvants within the same group:
  
  for i=1:length(expGroupPrefixes)
    currGroupNames = {};
    for j=1:length(expGroupStruct)
      currGroupNames = [currGroupNames  expGroupStruct(j).groupNames(i)];
    end
    %[figInd] = plotAllResponsesByGroup(mPeptResponses,groupNames,[currGroupNames, 'NC'],projectName,figInd,yLims,fontSize);
    [figInd] = plotSummaryResponsesByGroup(mPeptResponses,groupNames,[currGroupNames, 'NC'],projectName,summaryMethod,figInd,yLimsSummary,fontSize);

  end

  % plot post vs. pre results for all groups:
  [figInd] = plotPreVsPostSummaryResponsesByGroupObese(mPeptResponses,expGroupStruct,summaryMethod,figInd,yLimsSummary,fontSize);

end


% statistically compare and plot overall responses - summed across the array as boxplots of the four groups (pre/post WT/Ob):
% magnitude by ptid across entire array and the H7N9 array.
totalMagnitude = sum(mPeptResponses,2);
H7N9magnitude = sum(mPeptResponses(:,ShaInds),2);

for i=1:length(expGroupStruct)

  % compare groups using ranksum  test focusing only on post boost responses:

  % note that expGroupStruct is organized always in the same structure - this is why real indices (3,4) can safely be used here
  WT_postInds = expGroupStruct(i).groupInds{3};
  Ob_postInds = expGroupStruct(i).groupInds{4};

  % ranksum comparisons 
  [expGroupStruct(i).totalMagP] = ranksum(totalMagnitude(WT_postInds),totalMagnitude(Ob_postInds));
  [expGroupStruct(i).H7N9magP]  = ranksum(H7N9magnitude(WT_postInds),H7N9magnitude(Ob_postInds));

  % plot total magnitude by treatment group
  figure(figInd)
  a = gca;
  set(a,'FontSize',16)
  ptidInds   = expGroupStruct(i).inds;
  ptidLabels = expGroupStruct(i).labels(ptidInds);
  % trick to remove labels from empty groups for boxplot...,
  uniqLabels = unique(ptidLabels); 
  boxplot(H7N9magnitude(ptidInds),ptidLabels,'labels',treatmentLabels(uniqLabels))
  title(sprintf('%s H7N9 magnitude ranksum p = %f',expGroupStruct(i).adjuvant,expGroupStruct(i).H7N9magP));
  figInd = figInd + 1;
end



   %print figures to files:
if (printFlag)
    figPath = [savePath,filesep,'FIG',filesep];
    for j=1:figInd
        figName = [figPath,projectName,'Figure_',num2str(j)];
        [figHandle] = printFigure(j,figName);
        close(figure(j));
    end
end
%}

