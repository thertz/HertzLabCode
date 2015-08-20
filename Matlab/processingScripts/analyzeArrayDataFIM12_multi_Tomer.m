

% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

% To create new project files save as changing description, then search < > filling in with the appropriate info.
%
clear all;
close all;

ZachFlag  = 0;
clusterPlotFlag = 1;
simplePlotFlag  = 0;
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
%multiple dates can be entered, but must have a corresponding prefix
experimentDates = {'05_27_2014','06_02_2014','06_03_2014','06_04_2014','06_05_2014','06_13_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'FIM12_zm_05_27_14_','FIM12_zm_06_02_14_','FIM12_zm_06_03_14_','FIM12_zm_06_04_14_','FIM12_zm_06_05_14_','FIM12_zm_06_13_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

%experimentDates = {'05_27_2014','06_02_2014','06_03_2014','06_05_2014','06_13_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
%expPrefix    = {'FIM12_zm_05_27_14_','FIM12_zm_06_02_14_','FIM12_zm_06_03_14_','FIM12_zm_06_05_14_','FIM12_zm_06_13_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix


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
fontSize                 = 16;
%expGroups                = {'<3>','<2>','<1>'};% groups listed in sample to slide mapping
yLims                    = [0 65000];% max hight of y in graphs
minThreshold             = 10000; %minimal threshold for peak responses, used by findPeaks.
bgThreshold              = 5000; % threshold by which to flag antigens as ones that have high background and should
                                 % be removed from analysis (typically using a BSA negative control).
% Unsupervised clustering parameters:
numClusters          = 3; % number of clusters in the data (or expected number)
minResponseThreshold = 2000; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise

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

%transform all responses into log10 responses:
mPeptResponses(mPeptResponses<=0) = 100;
mPeptResponses = log10(mPeptResponses); 

antigenNames = arrayData(1).antigenNames;
strainNames  = {'Cal','Vic','Wis'};    
%index sets for the different strains/proteins:
Cal_HA_Inds = strmatch('Cal_HA',antigenNames);
Cal_NA_Inds = strmatch('Cal_NA',antigenNames);
Vic_HA_Inds = strmatch('Vic_HA',antigenNames);
Vic_NA_Inds = strmatch('Vic_NA',antigenNames);
Wis_HA_Inds = strmatch('Wis_HA',antigenNames);
Wis_NA_Inds = strmatch('Wis_NA',antigenNames);

CalInds = [Cal_HA_Inds;Cal_NA_Inds];
VicInds = [Vic_HA_Inds; Vic_NA_Inds];
WisInds = [Wis_HA_Inds; Wis_NA_Inds];

maskLabels = {'Cal_HA_Inds','Cal_NA_Inds','Vic_HA_Inds','Vic_NA_Inds','Wis_HA_Inds','Wis_NA_Inds'};
treatmentLabels = {'reg dose','high dose'};
infectionLabels = {'H3N2 neg','H3N2 pos'};

%load FIM12 labels!!!:
[FIM12struct] = parseFIM12dataFile(rootPath);

%remove negative control responses as measured by BSA alone:
BSAinds     = strmatch('BSA',ptids)';
nonBSAinds  = setdiff([1:length(ptids)],BSAinds);
bgResponses = max(mPeptResponses(BSAinds,:));
mPeptResponses(nonBSAinds,:) = mPeptResponses(nonBSAinds,:) - bgResponses;


%LIMIT ANALYSIS only to samples that are in the current ptid set - removing repeats but retaining BSA responses.
[C,IA,IB]      = intersect(ptids,{FIM12struct.slideID});
FIM12struct    = FIM12struct(IB);
mPeptResponses = mPeptResponses([IA,BSAinds],:);
ptids          = ptids([IA,BSAinds]);
ptidInds       = [1:length(IA)]; % INDICES NOT INCLUDING BSA responses

% write data to CSV file (patient samples only) for training classifiers.
rawDataMatrix = [mPeptResponses(1:length(IA),:) [FIM12struct.label]', [FIM12struct.H3N2pos]'];
rawDataFilename = [savePath,'rawDataFim12.csv'];
csvwrite(rawDataFilename, rawDataMatrix);

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

figInd = 1;

% magnitude by ptid
totalMagnitude = sum(mPeptResponses,2);
H3N2magnitude = sum(mPeptResponses(:,VicInds),2);

[totalMagP] = ranksum(totalMagnitude(regDoseInds),totalMagnitude(highDoseInds));
[H3N2magP] = ranksum(H3N2magnitude(regDoseInds),H3N2magnitude(highDoseInds));

[totalMagInfP] = ranksum(totalMagnitude(H3N2negInds),totalMagnitude(H3N2posInds));
[H3N2magInfP] = ranksum(H3N2magnitude(H3N2negInds),H3N2magnitude(H3N2posInds));

% plot total magnitude by treatment group
figure(1)
a = gca;
set(a,'FontSize',16)
boxplot(totalMagnitude(ptidInds),[FIM12struct.label],'labels',treatmentLabels)
title(sprintf('total magnitude ranksum p = %.3f',totalMagP));


figure(2)
a = gca;
set(a,'FontSize',16)
boxplot(H3N2magnitude(ptidInds),[FIM12struct.label],'labels',treatmentLabels);
title(sprintf('H3N2 magnitude ranksum p = %.3f',H3N2magP));


% plot total magnitude by infection status
figure(3)
a = gca;
set(a,'FontSize',16)
boxplot(totalMagnitude(ptidInds),[FIM12struct.H3N2pos],'labels',infectionLabels)
title(sprintf('total magnitude ranksum p = %.3f',totalMagInfP));


figure(4)
a = gca;
set(a,'FontSize',16)
boxplot(H3N2magnitude(ptidInds),[FIM12struct.H3N2pos],'labels',infectionLabels);
title(sprintf('H3N2 magnitude ranksum p = %.3f',H3N2magInfP));


figure(5)
a = gca;
set(a,'FontSize',16)
boxplot(H3N2magnitude(ptidInds),quadLabels,'labels',quadStrs);
title(sprintf('H3N2 magnitude ranksum p = %.3f',H3N2magInfP));

% breadth by ptid
posThreshold = log10(20000);
totalBreadth = sum(mPeptResponses>posThreshold,2);
H3N2breadth = sum(mPeptResponses(:,VicInds)>posThreshold,2);

[totalBreadthP] = ranksum(totalBreadth(regDoseInds),totalBreadth(highDoseInds));
[H3N2breadthP] = ranksum(H3N2breadth(regDoseInds),H3N2breadth(highDoseInds));

[totalbreadthInfP] = ranksum(totalBreadth(H3N2negInds),totalBreadth(H3N2posInds));
[H3N2breadthInfP] = ranksum(H3N2breadth(H3N2negInds),H3N2breadth(H3N2posInds));

% plot total breadth by treatment group
figure(6)
a = gca;
set(a,'FontSize',16)
boxplot(totalBreadth(ptidInds),[FIM12struct.label],'labels',treatmentLabels)
title(sprintf('total breadth ranksum p = %.3f',totalMagP));


figure(7)
a = gca;
set(a,'FontSize',16)
boxplot(H3N2breadth(ptidInds),[FIM12struct.label],'labels',treatmentLabels);
title(sprintf('H3N2 breadth ranksum p = %.3f',H3N2magP));


% plot total breadth by infection status
figure(8)
a = gca;
set(a,'FontSize',16)
boxplot(totalBreadth(ptidInds),[FIM12struct.H3N2pos],'labels',infectionLabels)
title(sprintf('total breadth ranksum p = %.3f',totalMagInfP));


figure(9)
a = gca;
set(a,'FontSize',16)
boxplot(H3N2breadth(ptidInds),[FIM12struct.H3N2pos],'labels',infectionLabels);
title(sprintf('H3N2 breadth ranksum p = %.3f',H3N2magInfP));


figure(10)
a = gca;
set(a,'FontSize',16)
boxplot(H3N2breadth(ptidInds),quadLabels,'labels',quadStrs);
title(sprintf('H3N2 breadth ranksum p = %.3f',H3N2magInfP));
keyboard;





figure(figInd);
for i=1:length(maskLabels)
  subplot(6,1,i)
  plotStr = ['bar(doseRespMat(',maskLabels{i},',:))'];
  eval(plotStr)
  legend('reg dose','high dose')
  title(strrep(maskLabels{i}(1:6),'_',' '));
end
figInd = figInd + 1;

% plot average responses by infection status
meanRespsH3N2pos = mean(mPeptResponses(H3N2posInds,:));
meanRespsH3N2neg = mean(mPeptResponses(H3N2negInds,:));
infectionRespMat = [meanRespsH3N2neg' meanRespsH3N2pos'];

figure(figInd);
for i=1:length(maskLabels)
  subplot(6,1,i)
  plotStr = ['bar(infectionRespMat(',maskLabels{i},',:))'];
  eval(plotStr)
  legend('infected','un-infected')
  title(strrep(maskLabels{i}(1:6),'_',' '));
end
figInd= figInd + 1 ;

%compare responses of all 4 gropus (infected/un-infected reg/high dose):
groupResp(1,:) = mean(mPeptResponses(aInds,:));
groupResp(2,:) = mean(mPeptResponses(bInds,:));
groupResp(3,:) = mean(mPeptResponses(cInds,:));
groupResp(4,:) = mean(mPeptResponses(dInds,:));
groupResp = groupResp'; % for plotting bar plots

figure(figInd);
for i=1:length(maskLabels)
  subplot(6,1,i)
  plotStr = ['bar(groupResp(',maskLabels{i},',:))'];
  eval(plotStr)
  legend('neg-reg','neg-high','pos-reg','pos-high');
  title(strrep(maskLabels{i}(1:6),'_',' '));
end
figInd= figInd + 1 ;

bar(groupResp')


if (clusterPlotFlag)   

%cluster by strain and/or protein
  i
%  for i=1:length(strainNames)
     %cluster data using unserupvised correlation clustering and plot cluster dendrogram (tree) - need to define the
     %following parameters up-top: numClusters and minThreshold
%    [distMat{i},Zstruct{i},clusterLabels{i}] = clusterSamplesByResponseVectors(mPeptResponses([,numClusters,minResponseThreshold);

  [distMatCal,ZstructCal,clusterLabelsCal] = clusterSamplesByResponseVectors(mPeptResponses(:,[Cal_HA_Inds ; Cal_NA_Inds]),numClusters,minResponseThreshold);
  [corrMatCal] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,[Cal_HA_Inds ; Cal_NA_Inds]),minResponseThreshold);

  [distMatVic,ZstructVic,clusterLabelsVic] = clusterSamplesByResponseVectors(mPeptResponses(:,[Vic_HA_Inds ; Vic_NA_Inds]),numClusters,minResponseThreshold);
  [corrMatVic] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,[Vic_HA_Inds ; Vic_NA_Inds]),minResponseThreshold);
n
  [distMatWis,ZstructWis,clusterLabelsWis] = clusterSamplesByResponseVectors(mPeptResponses(:,[Wis_HA_Inds ; Wis_NA_Inds]),numClusters,minResponseThreshold);
  [corrMatWis] = computeSampleCorrMatByResponseVectors(mPeptResponses(:,[Wis_HA_Inds ; Wis_NA_Inds]),minResponseThreshold);

  % cluster with a sepcific set of antigens:
  %[distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponses,numClusters,minResponseThreshold);
   
  % PLOT dendrogram:
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = 'Victoria H3N2'; %dendrogram title
  [figInd] = plotResponseDendrogram(ZstructVic,{FIM12struct.subtype},figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabelsVic);
  [figInd] = plotSampleCorrmat(figInd,corrMatVic(I,I),ptids(I),[0 1],fontSize);
  title(titleStr);

  titleStr = 'Cal H1N1'; %dendrogram title
  [figInd] = plotResponseDendrogram(ZstructCal,{FIM12struct.treatment},figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabelsCal);
  [figInd] = plotSampleCorrmat(figInd,corrMatCal(I,I),ptids(I),[0 1],fontSize);
  title(titleStr);

  titleStr = 'Wisconsin B'; %dendrogram title
  [figInd] = plotResponseDendrogram(ZstructWis,{FIM12struct.treatment},figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabelsWis);
  [figInd] = plotSampleCorrmat(figInd,corrMatWis(I,I),ptids(I),[0 1],fontSize);
  title(titleStr);
end  

%{
-getting error from this chunk- does it need to be flagged at the top?
figure(figInd);
for i=1:numClusters


  gInds = find(clusterLabelsVic == i);
  mMat(i,:) = mean(mPeptResponses(gInds,[Vic_HA_Inds ; Vic_NA_Inds]));

  subplot(numClusters,1,i);
  plot(mPeptResponses(gInds,[Vic_HA_Inds ; Vic_NA_Inds])');
  title('Victoria H3N2');

end

figure;
bar(mMat(:,[1:length(Vic_HA_Inds)])');
title('Vic HA');

figure;
bar(mMat(:,[length(Vic_HA_Inds)+1:end])');
title('Vic NA');
%}

  %-------------------------------------------------------------------------------%  
  % plotting functions
  %-------------------------------------------------------------------------------%
if (simplePlotFlag)
  [figInd] = plotAllResponsesByGroup(mPeptResponses,groupNames,expGroups,projectName,figInd,yLims,fontSize);%
  [figInd] = plotAllResponsesByPtid(mPeptResponses,ptids,'<>',projectName,figInd,yLims)

end


   %print figures to files:
if (printFlag)
    figPath = [savePath,filesep,'FIG',filesep];
    for j=1:figInd
        figName = [figPath,projectName,'Figure_',num2str(j)];
        printFig(figure(j),figName,30,24);
    end
end


