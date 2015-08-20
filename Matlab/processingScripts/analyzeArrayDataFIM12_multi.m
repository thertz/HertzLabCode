

% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

%
clear all;
close all;

ZachFlag  = 0;
plotFlag  = 1;
printFlag = 0;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'FIM12';

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
experimentDates = {'05_27_2014','06_02_2014','06_03_2014','06_04_2014','06_05_2014','06_13_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'FIM12_zm_05_27_14_','FIM12_zm_06_02_14_','FIM12_zm_06_03_14_','FIM12_zm_06_04_14_','FIM12_zm_06_05_14_','FIM12_zm_06_13_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix



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
minThreshold = 2048 %minimal threshold for peak responses, used by findPeaks.


% Unsupervised clustering parameters:
numClusters          = 2; % number of clusters in the data (or expected number)
minResponseThreshold = 2048; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise
bgThreshold          = 4096; % threshold by which to flag antigens as ones that have high background and should
                             % be removed from analysis (typically using a BSA negative control).

%-------------------------------------------------------------------------------%
% Stat comparison parameters: (also project specific)
%-------------------------------------------------------------------------------%
% Filters used for treatment-blinded filtering of responses to single antigens used for statistical analysis.
% currently only one filter is implemented which is % responders within the entire treatment set (ignoring all controls) 
filterNames = {'PercentResponders'};
filterThresholds = {0.2};
treatmentGroups = {'reg-dose','high-dose'}; % Names of the two groups that are to be compared statistically below (currently supports only pairwise comparisons)

H3N2_antigenic7_SiteInds =  [160, 170, 171, 173, 174, 204, 208];

% data Labels
maskLabels = {'Cal_HA_Inds','Cal_NA_Inds','Vic_HA_Inds','Vic_NA_Inds','Wis_HA_Inds','Wis_NA_Inds'};
treatmentLabels = {'reg dose','high dose'};
infectionLabels = {'H3N2 neg','H3N2 pos'};

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
%remove negative control responses as measured by BSA alone:
BSAinds     = strmatch('BSA',ptids)';
nonBSAinds  = setdiff([1:length(ptids)],BSAinds);
bgResponses = max(mPeptResponses(BSAinds,:));

% split out BSA responses - to allow label consistency with the FIM12struct
BSAresponses   = mPeptResponses(BSAinds,:);

slideNames = slideNames(nonBSAinds);
ptids = ptids(nonBSAinds);
groupNames = groupNames(nonBSAinds);

mPeptResponses = mPeptResponses(nonBSAinds,:);
mPeptResponses = mPeptResponses - repmat(bgResponses,length(nonBSAinds),1);

% clean data and take log2 transform:
mPeptResponses(mPeptResponses<1024) = 1024;
%mPeptResponses = log2(mPeptResponses);


antigenNames = arrayData(1).antigenNames;
strainNames  = {'Cal','Vic','Wis'};    
%index sets for the different strains/proteins:
Cal_HA_Inds = strmatch('Cal_HA',antigenNames);
Cal_NA_Inds = strmatch('Cal_NA',antigenNames);
Vic_HA_Inds = strmatch('Vic_HA',antigenNames);
Vic_NA_Inds = strmatch('Vic_NA',antigenNames);
Wis_HA_Inds = strmatch('Wis_HA',antigenNames);
Wis_NA_Inds = strmatch('Wis_NA',antigenNames);

% Struct use to encode all antigen masks to allwo streamlining code for iterating over mask sets
antigenMasksStruct(1).name = maskLabels{1};
antigenMasksStruct(1).inds = Cal_HA_Inds;

antigenMasksStruct(2).name = maskLabels{2};
antigenMasksStruct(2).inds = Cal_NA_Inds;

antigenMasksStruct(3).name = maskLabels{3};
antigenMasksStruct(3).inds = Vic_HA_Inds;

antigenMasksStruct(4).name = maskLabels{4};
antigenMasksStruct(4).inds = Vic_NA_Inds;

antigenMasksStruct(5).name = maskLabels{5};
antigenMasksStruct(5).inds = Wis_HA_Inds;

antigenMasksStruct(6).name = maskLabels{6};
antigenMasksStruct(6).inds = Wis_NA_Inds;

CalInds = [Cal_HA_Inds;Cal_NA_Inds];
VicInds = [Vic_HA_Inds; Vic_NA_Inds];
WisInds = [Wis_HA_Inds; Wis_NA_Inds];


%load FIM12 labels!!!:
% Based on case-control dataset sent by Andrew Dunning
[FIM12struct] = parseFIM12caseControlDataFile(rootPath)


% order FIM12 struct by current ptids, and remove ones that were not processed:
ptidInds = [];
remInds  = []; % indices of ptids to be removed since they do not appear in case-control set sent by Dunning
for i=1:length(ptids)
  currInd = strmatch(ptids{i},{FIM12struct.slideID});
  
  % this indicates major issue with sample labeling 
  if(isempty(currInd))
    if(~strcmp(ptids{i},'BSA'))
      disp(sprintf('Error: cannot find label for ptid %s',ptids{i}))
      remInds = [remInds i];
    end
  end

  ptidInds = [ptidInds currInd];
end

%remove all ptids that were not in the case-control set (why are there some?)
slideNames(remInds) = [];
ptids(remInds) = [];
groupNames(remInds) = [];
mPeptResponses(remInds,:) = [];

% reorder FIM12struct by current ptids
FIM12struct = FIM12struct(ptidInds);


% load FIM12 peptide data:
[pepData] = readFLU_FIM12_PeptideData(rootPath);

% ALTERNATE VERSION BASED ON FULL LABEL FILE - % WARNING - CURRENTLY BROKEN DUE TO REMOVING BSA SAMPLES FROM DATA
%[FIM12struct] = parseFIM12dataFile(rootPath);

% Limit analyses only to samples that are in the current ptid set - removing repeats but retaining BSA responses. - for alternate version of FIM12struct - 

% [C,IA,IB]      = intersect(ptids,{FIM12struct.slideID});
% FIM12struct    = FIM12struct(IB);
% mPeptResponses = mPeptResponses([IA,BSAinds],:);
% ptids          = ptids([IA,BSAinds]);
% ptidInds       = [1:length(IA)]; % INDICES NOT INCLUDING BSA responses

% write data to CSV file (patient samples only) for training classifiers.
rawDataMatrix = [mPeptResponses [FIM12struct.label]', [FIM12struct.H3N2pos]'];
rawDataFilename = [savePath,'rawDataFim12.csv'];
%csvwrite(rawDataFilename, rawDataMatrix);

%-------------------------------------------------------------------------------%
% analysis 
%-------------------------------------------------------------------------------%
%find peak responses to visually verify these on images:
[peakResps] = findPeakResponses(mPeptResponses,minThreshold);

regDoseInds  = find([FIM12struct.label]==0);
highDoseInds = find([FIM12struct.label]==1);

H3N2posInds  = find([FIM12struct.H3N2pos]);
H3N2negInds  = find(~[FIM12struct.H3N2pos]);

aInds = intersect(H3N2negInds,regDoseInds);
bInds = intersect(H3N2negInds,highDoseInds);
cInds = intersect(H3N2posInds,regDoseInds);
dInds = intersect(H3N2posInds,highDoseInds);

quadLabels = zeros(1,length(ptidInds));
quadLabels(aInds) = 1;
quadLabels(bInds) = 2;
quadLabels(cInds) = 3;
quadLabels(dInds) = 4;
quadStrs = {'neg-reg','neg-high','pos-reg','pos-high'};


%-------------------------------------------------------------------------------%
% Compute pairwise statistics comparing treatment groups of interest
% 
%  
%-------------------------------------------------------------------------------%

% Reg-dose vs. High-dose:
respMat = mPeptResponses;%([regDoseInds,highDoseInds],:);
treatmentLabelsVec  = [FIM12struct.label];

[antigenStatsTreatment] = compareGroupsByAntigenUsingAntigenMasks(respMat,treatmentLabelsVec,antigenMasksStruct,filterNames,filterThresholds,minResponseThreshold);

% test to see if the H3 sites are in the antigenic 7 filter
antigenic7_Inds = [];
H3inds = strmatch('Vic_HA',{pepData.name});
for i = 1:length(H3inds)
    if(~isempty(intersect(H3N2_antigenic7_SiteInds,[pepData(H3inds(i)).begInd:pepData(H3inds(i)).endInd])))
      antigenic7_Inds = [antigenic7_Inds, i];
    end
end
antigenic7_pValues = antigenStatsTreatment(3).singleAntigen_pValuesRankSum(antigenic7_Inds);
% Q-value correction - adjust all p-values that are within the set of filtered inds
[FDR, qValues] = mafdr(antigenic7_pValues);


[antigenStatsTreatment(3)] =  ;


% OUTPUT results to command line:
disp(sprintf('Results for comparisons between treatment groups (reg-dose vs. high-dose'));
disp(sprintf('*----------------------------------------------------------------------------*'))
for i=1:length(antigenStatsTreatment)
  disp(sprintf('Comparing overall responses across each entire vaccine antigen %s', antigenStatsTreatment(i).name))
  disp(sprintf('Ranksum p-value: %0.5g\n\n',antigenStatsTreatment(i).totResponses_pValue))
end
disp(sprintf('*----------------------------------------------------------------------------*'))

disp(sprintf('Results for comparisons of single antigens between treatment groups (reg-dose vs. high-dose'));
disp(sprintf('*----------------------------------------------------------------------------*'))
for i=1:length(antigenStatsTreatment)
  disp(sprintf('Comparing single antigen responses for %s:',antigenStatsTreatment(i).name ))
  for j=1:length(antigenStatsTreatment(i).qIndsRankSum)
    disp(sprintf('q-value: %0.5g\t antigen # %d',antigenStatsTreatment(i).qValuesRankSum(j),antigenStatsTreatment(i).qIndsRankSum(j)));
  end
  disp(sprintf('\n'))
end
disp(sprintf('*----------------------------------------------------------------------------*'))



% Infected vs. un-infected
respMat = mPeptResponses;%([H3N2posInds,H3N2negInds],:);
infectionLabelsVec  = [FIM12struct.H3N2pos];
[antigenStatsInfection] = compareGroupsByAntigenUsingAntigenMasks(respMat,infectionLabelsVec,antigenMasksStruct,filterNames,filterThresholds,minResponseThreshold);

% OUTPUT results to command line:
disp(sprintf('Results for comparisons between outcome groups (infected vs. un-infected'));
disp(sprintf('*----------------------------------------------------------------------------*'))
for i=1:length(antigenStatsInfection)
  disp(sprintf('Comparing overall responses across each entire vaccine antigen %s', antigenStatsInfection(i).name))
  disp(sprintf('Ranksum p-value: %0.5g\n\n',antigenStatsInfection(i).totResponses_pValue))
end
disp(sprintf('*----------------------------------------------------------------------------*'))

disp(sprintf('Results for comparisons of single antigens between treatment groups (reg-dose vs. high-dose'));
disp(sprintf('*----------------------------------------------------------------------------*'))
for i=1:length(antigenStatsInfection)
  disp(sprintf('Comparing single antigen responses for %s:',antigenStatsTreatment(i).name ))
  for j=1:length(antigenStatsInfection(i).qIndsRankSum)
    disp(sprintf('q-value: %0.5g\t antigen # %d',antigenStatsInfection(i).qValuesRankSum(j),antigenStatsInfection(i).qIndsRankSum(j)));
  end
  disp(sprintf('\n'))
end
disp(sprintf('*----------------------------------------------------------------------------*'))



%-------------------------------------------------------------------------------
% magnitude by ptid
%
% Compare total magnitude of responses by treatment assignment and infection status.
%-------------------------------------------------------------------------------

totalMagnitude = sum(mPeptResponses,2);
H3N2magnitude = sum(mPeptResponses(:,VicInds),2);

[totalMagnitudeTreatmentP] = ranksum(totalMagnitude(regDoseInds),totalMagnitude(highDoseInds));
[totalMagnitudeInfP] = ranksum(totalMagnitude(H3N2negInds),totalMagnitude(H3N2posInds));

[H3N2_MagnitudeTreatmentP] = ranksum(H3N2magnitude(regDoseInds),H3N2magnitude(highDoseInds));
[H3N2_MagnitudeInfP] = ranksum(H3N2magnitude(H3N2negInds),H3N2magnitude(H3N2posInds));

statStruct.totalMagnitudeTreatmentP = totalMagnitudeTreatmentP;
statStruct.totalMagnitudeInfP       = totalMagnitudeInfP;

statStruct.H3N2_MagnitudeTreatmentP = H3N2_MagnitudeTreatmentP;
statStruct.H3N2_MagnitudeInfP       = H3N2_MagnitudeInfP;



% plot total magnitude by treatment group
if(plotFlag)
  figure(1);
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(totalMagnitude,[FIM12struct.label],[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',treatmentLabels)
  title(sprintf('total magnitude ranksum p = %.3f',totalMagnitudeTreatmentP));


  figure(2)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(H3N2magnitude,[FIM12struct.label],[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',treatmentLabels)
  title(sprintf('H3N2 magnitude ranksum p = %.3f',H3N2_MagnitudeTreatmentP));


  % plot total magnitude by infection status
  figure(3)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(totalMagnitude,[FIM12struct.H3N2pos],[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',infectionLabels)
  title(sprintf('total magnitude ranksum p = %.3f',totalMagnitudeInfP));


  figure(4)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(H3N2magnitude,[FIM12struct.H3N2pos],[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',infectionLabels)
  title(sprintf('H3N2 magnitude ranksum p = %.3f',H3N2_MagnitudeInfP));
  

  figure(5)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(H3N2magnitude,quadLabels,[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',quadStrs)
  title(sprintf('H3N2 magnitude by groups'));
end


%-------------------------------------------------------------------------------
% breadth by ptid
%
% define arbitrary low threshold for positive response and compare the breadth 
% of responses of groups by treatment assignemnt and infection status.
%-------------------------------------------------------------------------------

%posThreshold = log2(2048);
posThreshold = 4096;

totalBreadth = sum(mPeptResponses>posThreshold,2);
H3N2_Breadth = sum(mPeptResponses(:,VicInds)>posThreshold,2);

[totalBreadthTreatmentP] = ranksum(totalBreadth(regDoseInds),totalBreadth(highDoseInds));
[totalbreadthInfP] = ranksum(totalBreadth(H3N2negInds),totalBreadth(H3N2posInds));

[H3N2_BreadthTreatmentP]  = ranksum(H3N2_Breadth(regDoseInds),H3N2_Breadth(highDoseInds));
[H3N2_BreadthInfP]  = ranksum(H3N2_Breadth(H3N2negInds),H3N2_Breadth(H3N2posInds));

statStruct.totalBreadthTreatmentP = totalBreadthTreatmentP;
statStruct.totalbreadthInfP       = totalbreadthInfP;

statStruct.H3N2_BreadthTreatmentP = H3N2_BreadthTreatmentP;
statStruct.H3N2_BreadthInfP       = H3N2_BreadthInfP;


% plot total breadth by treatment group
if(plotFlag)
  figure(6)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(totalBreadth,[FIM12struct.label],[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',treatmentGroups);
  title(sprintf('total breadth ranksum p = %.3f',totalBreadthTreatmentP));


  figure(7)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(H3N2_Breadth,[FIM12struct.label],[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',treatmentGroups);
  title(sprintf('H3N2 breadth ranksum p = %.3f',H3N2_BreadthTreatmentP));

  % plot total breadth by infection status
  figure(8)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(totalBreadth,[FIM12struct.H3N2pos],[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',infectionLabels);
  title(sprintf('total breadth ranksum p = %.3f',totalbreadthInfP));



  figure(9)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(H3N2_Breadth,[FIM12struct.H3N2pos],[],'sdline')
  set(a,'FontSize',fontSize)
  set(a,'XtickLabel',infectionLabels);
  title(sprintf('H3N2 breadth ranksum p = %.3f',H3N2_BreadthInfP));


  figure(10)
  a = gca;
  set(a,'FontSize',fontSize)
  h = notBoxPlot(H3N2_Breadth,quadLabels,[],'sdline')
  set(a,'FontSize',fontSize);
  set(a,'XtickLabel',quadStrs);
  title(sprintf('H3N2 breadth by groups'));
  
end

%-------------------------------------------------------------
% analyze and plot average responses by treatment label
%-------------------------------------------------------------
meanRespsRegDose  = mean(mPeptResponses(regDoseInds,:));
meanRespsHighDose = mean(mPeptResponses(highDoseInds,:));
doseRespMat  = [meanRespsRegDose' meanRespsHighDose'];

% paired t-test on responses to each antigen by treatment assignment
[H,P] = ttest(meanRespsRegDose,meanRespsHighDose);
[H1,P1] = ttest(meanRespsRegDose(VicInds),meanRespsHighDose(VicInds));

statStruct.meanRespByTreatmentP = P;
statStruct.meanRespByTreatment_H3N2_P = P1;

figInd = 11;

figure(figInd);
for i=1:length(maskLabels)
  subplot(6,1,i)
  plotStr = ['bar(doseRespMat(',maskLabels{i},',:))'];
  eval(plotStr)
  legend('reg dose','high dose')
  title(strrep(maskLabels{i}(1:6),'_',' '));
end
figInd = figInd + 1;

%-------------------------------------------------------------
% analyze and plot average responses by infection status
%-------------------------------------------------------------
meanRespsH3N2pos = mean(mPeptResponses(H3N2posInds,:));
meanRespsH3N2neg = mean(mPeptResponses(H3N2negInds,:));
infectionRespMat = [meanRespsH3N2pos' meanRespsH3N2neg'];

% paired t-test on responses to each antigen by infection status
[H,P] = ttest(meanRespsH3N2pos,meanRespsH3N2neg);
[H1,P1] = ttest(meanRespsH3N2pos(VicInds),meanRespsH3N2neg(VicInds));

statStruct.meanRespByInfectionStatusP = P;
statStruct.meanRespByInfectionStatus_H3N2_P = P1;


figure(figInd);
for i=1:length(maskLabels)
  subplot(6,1,i)
  plotStr = ['bar(infectionRespMat(',maskLabels{i},',:))'];
  eval(plotStr)
  legend('infected','un-infected')
  title(strrep(maskLabels{i}(1:6),'_',' '));
end
figInd= figInd + 1 ;


%-------------------------------------------------------------------------
%compare responses of all 4 gropus (infected/un-infected reg/high dose):
%-------------------------------------------------------------------------
groupResp(1,:) = mean(mPeptResponses(aInds,:)); % neg-reg
groupResp(2,:) = mean(mPeptResponses(bInds,:)); % neg-high
groupResp(3,:) = mean(mPeptResponses(cInds,:)); % pos-reg 
groupResp(4,:) = mean(mPeptResponses(dInds,:)); % pos-high
groupResp = groupResp'; % for plotting bar plots


% compare groups using ranksum test (may require stat correction for structure)
indLabels = {'aInds','bInds','cInds','dInds'};
for i=1:4
  for j=i+1:4
    cmdStr= ['[H3N2_P(i,j), H3N2_H(i,j)] = ranksum(H3N2magnitude(',indLabels{i},',:),H3N2magnitude(',indLabels{j},',:))'];
    eval(cmdStr);

    cmdStr1= ['[P(i,j), H(i,j)] = ranksum(totalMagnitude(',indLabels{i},',:),totalMagnitude(',indLabels{j},',:))'];
    eval(cmdStr1);
  end
end


statStruct.groupRespP = P;
statStruct.gropuResp_H3N3_P = H3N2_P;


if(plotFlag)
  figure(figInd);
  for i=1:length(maskLabels)
    subplot(6,1,i)
    plotStr = ['bar(groupResp(',maskLabels{i},',:))'];
    eval(plotStr)
    legend('neg-reg','neg-high','pos-reg','pos-high');

    title(strrep(maskLabels{i}(1:6),'_',' '));
  end
  figInd= figInd + 1 ;
end

%-------------------------------------------------------------------------
% Unsupervised Clustering analyses
%-------------------------------------------------------------------------

  
%cluster by strain and/or protein
%following parameters up-top: numClusters and minThreshold
  
[distMatCal,ZstructCal,clusterLabelsCal] = clusterSamplesByResponseVectors(mPeptResponses(:,[Cal_HA_Inds ; Cal_NA_Inds]),numClusters,minResponseThreshold);
[corrMatCal] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,[Cal_HA_Inds ; Cal_NA_Inds]),minResponseThreshold);

[distMatVic,ZstructVic,clusterLabelsVic] = clusterSamplesByResponseVectors(mPeptResponses(:,[Vic_HA_Inds ; Vic_NA_Inds]),numClusters,minResponseThreshold);
[corrMatVic] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,[Vic_HA_Inds ; Vic_NA_Inds]),minResponseThreshold);

[distMatWis,ZstructWis,clusterLabelsWis] = clusterSamplesByResponseVectors(mPeptResponses(:,[Wis_HA_Inds ; Wis_NA_Inds]),numClusters,minResponseThreshold);
[corrMatWis] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,[Wis_HA_Inds ; Wis_NA_Inds]),minResponseThreshold);


% analyze clusters - compute Fisher's test within each cluster and also compute average responses within each cluster.
for i=1:numClusters
  inds = find(clusterLabelsVic == i);

  % compute Fisher's exact test on the 4 quadrants in each cluster. 
  a = length(intersect(find([FIM12struct(inds).H3N2pos]==0),find([FIM12struct(inds).label]==0)));
  b = length(intersect(find([FIM12struct(inds).H3N2pos]==0),find([FIM12struct(inds).label]==1)));
  c = length(intersect(find([FIM12struct(inds).H3N2pos]==1),find([FIM12struct(inds).label]==0)));
  d = length(intersect(find([FIM12struct(inds).H3N2pos]==1),find([FIM12struct(inds).label]==1)));
  [Ppos,Pneg,Pboth(i)]=Fisherextest(a,b,c,d);

  l(i) = length(inds); % used for deciding which clusters to consider - only consider clusters with 5 or more datapoints

  
  s(i)= sum([FIM12struct(inds).H3N2pos])./length(inds); % percent of H3N2 cases in cluster
  s1(i)= sum([FIM12struct(inds).label])./length(inds);  % percent of high-dose treatment in cluster
  
  mMat(i,:) = mean(mPeptResponses(inds,:)); % average response of members within a given cluster
end


if(plotFlag)
  % now plot the average response of each cluster of size greater than 5 to Vic antigens
  clusterInds = find(l>5);
  currL = length(clusterInds);
  for i=1:currL
    figure(figInd);
    subplot(currL,1,i);
    currInd = clusterInds(i);
    bar(mMat(currInd,Vic_HA_Inds)')
    a = gca;
    set(a,'Ylim',[0 40000])
    title(sprintf('Cluster #%d HA Vic Responses (n=%d), reg-dose : %.2f  ',currInd,l(currInd),s1(currInd)));

    figure(figInd+1)
    subplot(currL,1,i);
    bar(mMat(clusterInds(i),Vic_NA_Inds)')
    a = gca;
    set(a,'Ylim',[0 40000])
    title(sprintf('Cluster #%d NA Vic Responses (n=%d), reg-dose: %.2f  ',currInd,l(currInd),s1(currInd)));

  end
  figInd = figInd +2 ;

  clusterFontSize = 8;
  
  % PLOT dendrogram:
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = 'Victoria H3N2'; %dendrogram title
  [figInd] = plotResponseDendrogram(ZstructVic,{FIM12struct.treatment},figInd,clusterFontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabelsVic);
  [figInd] = plotSampleCorrmat(figInd,corrMatVic(I,I),{FIM12struct(I).treatment},[0 1],clusterFontSize);
  title(titleStr);

  % Plot single antigen response histograms:

% responses by treatment group: reg-dose vs. high-dose
  for i=1:length(antigenStatsTreatment)
    currInds = antigenMasksStruct(i).inds(antigenStatsTreatment(i).qIndsRankSum);
    for j=1:length(currInds)
      [figInd] = plotIndividualResponsesForSingleAntigensByGroup(respMat(:,currInds(j)),treatmentLabelsVec,treatmentGroups,minResponseThreshold,figInd,fontSize)
      [figInd] = plotSingleAntigenCDFsByGroup(respMat(:,currInds(j)),treatmentLabelsVec,treatmentGroups,minResponseThreshold,figInd,fontSize)

      title(sprintf('CDF for peptide %s , start-pos: %d, seq: %s, q-value = %.3f',strrep(pepData(currInds(j)).name,'_','-'),...
        pepData(currInds(j)).begInd,pepData(currInds(j)).sequence,antigenStatsTreatment(i).qValuesRankSum(j)));  
    end
  end

for i=1:length(antigenStatsInfection)
  currInds = antigenMasksStruct(i).inds(antigenStatsInfection(i).qIndsRankSum);
  for j=1:length(currInds)
    [figInd] = plotIndividualResponsesForSingleAntigensByGroup(respMat(:,currInds(j)),treatmentLabelsVec,treatmentGroups,minResponseThreshold,figInd,fontSize)
  end
end

  % titleStr = 'Cal H1N1'; %dendrogram title
  % [figInd] = plotResponseDendrogram(ZstructCal,{FIM12struct(I).case_control},figInd,clusterFontSize,bwFlag,titleStr);
  
  % %plot correlation matrix:    
  % [Y I] = sort(clusterLabelsCal);
  % [figInd] = plotSampleCorrmat(figInd,corrMatCal(I,I),{FIM12struct(I).case_control},[0 1],clusterFontSize);
  % title(titleStr);

  % titleStr = 'Wisconsin B'; %dendrogram title
  % [figInd] = plotResponseDendrogram(ZstructWis,{FIM12struct(I).case_control},figInd,clusterFontSize,bwFlag,titleStr);
  
  % %plot correlation matrix:    
  % [Y I] = sort(clusterLabelsWis);
  % [figInd] = plotSampleCorrmat(figInd,corrMatWis(I,I),{FIM12struct(I).case_control},[0 1],clusterFontSize);
  % title(titleStr);
end  

% plot raw data:
subplot(2,1,1)
plot(mPeptResponses(:,Vic_HA_Inds)')
a = gca;
set(a,'FontSize',fontSize);
title('Vic HA - FIM12 samples')
subplot(2,1,2)
plot(BSAresponses(:,Vic_HA_Inds)');
a = gca;
set(a,'FontSize',fontSize);
title('Vic HA - Negative Controls  (BSA)');

% plot raw data:
figInd = figInd + 1;
figure(figInd);
subplot(2,1,1)
bar(mean(mPeptResponses(:,Vic_HA_Inds)))
a = gca;
set(a,'FontSize',fontSize);
title('Vic HA average response- FIM12 samples')
subplot(2,1,2)
bar(mean(BSAresponses(:,Vic_HA_Inds)));
a = gca;
set(a,'FontSize',fontSize);
title('Vic HA avereage response - Negative Controls  (BSA)');



   %print figures to files:
if (printFlag)
    figPath = [savePath,'FIG',filesep];
    for j=1:10%figInd
        figName = [figPath,projectName,'_Figure_',num2str(j)]
        printFig(figure(j),figName,30,24);
        close(figure(j));
    end
end


