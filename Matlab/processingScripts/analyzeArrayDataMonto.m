

% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

% To create new project files save as changing description, then search < > filling in with the appropriate info.
% Monto data from U. Michigan, pre,post,end-of-season samples of individuals vaccinated with either flu mist or inactivated vaccine.
% SAMPLES were run on the FIM12 slides.
clear all;
close all;

ZachFlag  = 0;
plotFlag  = 1;
printFlag = 1;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'Monto';

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
experimentDates = {'07_15_2014','07_16_2014','07_23_2014','08_03_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'FIM12_zm_07_15_14_','FIM12_zm_07_16_14_','FIM12_zm_07_23_14_','FIM12_zm_08_03_14_'}; % all slides and slideToSampleMapping file for each date have the exact same filePrefix

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
fontSize                 = 10;
%expGroups                = {'<3>','<2>','<1>'};% groups listed in sample to slide mapping
yLims                    = [0 65000];% max hight of y in graphs
yLimsSummary             = [0 32500]; %max height of y in summary stat graphs
minThreshold             = 10000; %minimal threshold for peak responses, used by findPeaks.
bgThreshold              = 5000; % threshold by which to flag antigens as ones that have high background and should
                                 % be removed from analysis (typically using a BSA negative control).
% Unsupervised clustering parameters:
numClusters          = 6; % number of clusters in the data (or expected number)
minResponseThreshold = 2000; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise
summaryMethod = 'median';  % summary stat used for group comparisons and for plotting gropu responses - and cluster responses
%-------------------------------------------------------------------------------%
%defining elements outside of loop to add loop data too
slideNames      ={};
ptids           ={};
groupNames      ={};
mPeptResponses  =[];
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

%transform all responses into log10 responses:
%mPeptResponses(mPeptResponses<=0) = 100;
%mPeptResponses = log10(mPeptResponses); 

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
treatmentLabels = {'Inactivated_Fluzone','Live_Flumist'};
timepointPrefixes = {'Vac_','Post_','LT_'};

%remove negative control responses as measured by BSA alone:
BSAinds     = strmatch('BSA',ptids)';
nonBSAinds  = setdiff([1:length(ptids)],BSAinds);
bgResponses = median(mPeptResponses(BSAinds,:));
mPeptResponses(nonBSAinds,:) = mPeptResponses(nonBSAinds,:) - repmat(bgResponses,length(nonBSAinds),1);

% load peptide data:
[pepData] = readFLU_FIM12_PeptideData(rootPath);

% load ptid information including HI titers.
[ptidStruct] = parseMontoData(rootPath);

% now compute post and LT responses subtracting out BASELINE responses for each ptid separately!!!!
remInds      = []  % remove ptids for which not post or lt responses are available.
for i=1:length(ptidStruct)

  preInd  = strmatch([ptidStruct(i).ptid,'_1'],ptids);
  postInd = strmatch([ptidStruct(i).ptid,'_2'],ptids);
  ltInd   = strmatch([ptidStruct(i).ptid,'_3'],ptids);

  if( isempty(preInd) || isempty(postInd) || isempty(ltInd) )
    remInds = [remInds i];
    continue;
  else

    ptidStruct(i).preInd = preInd
    ptidStruct(i).postInd = postInd
    ptidStruct(i).ltInd = ltInd;
    ptidStruct(i).treatment = groupNames{preInd}(5:end); % remove Vac_ prefix.
    
    % sanity check - see that treatment as coded in our files matches the one in the original data manifest
    if(~strcmp(strrep(ptidStruct(i).treatmentOrig,' - ','_'),ptidStruct(i).treatment))
      disp('error in ptid treatment coding');
    end

    ptidStruct(i).postResponse = mPeptResponses(postInd,:);
    ptidStruct(i).ltResponse   = mPeptResponses(ltInd,:);

    ptidRespMatPre(i,:)  = mPeptResponses(preInd,:);
    ptidRespMatPost(i,:) = ptidStruct(i).postResponse - mPeptResponses(preInd,:);
    ptidRespMatLT(i,:)   = ptidStruct(i).ltResponse   - mPeptResponses(preInd,:);
  end
end

ptidStruct(remInds)        = [];
ptidRespMatPre(remInds,:)  = [];
ptidRespMatPost(remInds,:) = [];
ptidRespMatLT(remInds,:)   = [];

ptidLabels                 = {ptidStruct.treatment};
FlumistInds = strmatch('Live_Flumist',ptidLabels);
FluzoneInds = strmatch('Inactivated_Fluzone',ptidLabels);
labels(FluzoneInds) = 0;
labels(FlumistInds) = 1;

%-------------------------------------------------------------------------------%
% analysis 
%-------------------------------------------------------------------------------%

%-------------------------------------------------------------------------------%
% Unsupervised clustering
% 
%  Cluster data using unserupvised correlation clustering and plot cluster dendrogram (tree) - need to define the
%  following parameters up-top: numClusters and minThreshold
%-------------------------------------------------------------------------------%

% Cluster based on responses to the H3 - including HA and NA:    
[distMatVicPre,ZstructVicPre,clusterLabelsVicPre] = clusterSamplesByResponseVectors(ptidRespMatPre(:,[Vic_HA_Inds ; Vic_NA_Inds]),numClusters,minResponseThreshold);
[distMatVicPost,ZstructVicPost,clusterLabelsVicPost] = clusterSamplesByResponseVectors(ptidRespMatPost(:,[Vic_HA_Inds ; Vic_NA_Inds]),numClusters,minResponseThreshold);

% compute correlation matrices based on responses to H3:
[corrMatVicPre] = computeSampleCorrMatByResponseVectors(ptidRespMatPre(:,[Vic_HA_Inds ; Vic_NA_Inds]),minResponseThreshold);
[corrMatVicPost] = computeSampleCorrMatByResponseVectors(ptidRespMatPost(:,[Vic_HA_Inds ; Vic_NA_Inds]),minResponseThreshold);


%-------------------------------------------------------------------------------%
% Statistical group comparisons
%-------------------------------------------------------------------------------%

% magnitude by ptid - pre/post/lt
totalMagnitudePre = sum(ptidRespMatPre,2);
H3N2magnitudePre  = sum(ptidRespMatPre(:,VicInds),2);

[totalMagPreP] = ranksum(totalMagnitudePre(FlumistInds),totalMagnitudePre(FluzoneInds));
[H3N2magPreP] = ranksum(H3N2magnitudePre(FlumistInds),H3N2magnitudePre(FluzoneInds));

totalMagnitudePost = sum(ptidRespMatPost,2);
H3N2magnitudePost  = sum(ptidRespMatPost(:,VicInds),2);

[totalMagPostP] = ranksum(totalMagnitudePost(FlumistInds),totalMagnitudePost(FluzoneInds));
[H3N2magPostP] = ranksum(H3N2magnitudePost(FlumistInds),H3N2magnitudePost(FluzoneInds));


[r,p] = corr(log2(double([ptidStruct.H3_Wis_HI_preTiter]))',H3N2magnitudePre,'type','Spearman');
[r,p] = corr(log2(double([ptidStruct.H3_Wis_HI_postTiter]))',H3N2magnitudePre,'type','Spearman');


% breadth by ptid
posThreshold = 10000;
totalBreadthPre = sum(ptidRespMatPre>posThreshold,2);
H3N2breadthPre = sum(ptidRespMatPre(:,VicInds)>posThreshold,2);
  
totalBreadthPost = sum(ptidRespMatPost>posThreshold,2);
H3N2breadthPost = sum(ptidRespMatPost(:,VicInds)>posThreshold,2);
   
[totalBreadthPreP] = ranksum(totalBreadthPre(FluzoneInds),totalBreadthPre(FlumistInds));
[H3N2breadthPreP]  = ranksum(H3N2breadthPre(FluzoneInds),H3N2breadthPre(FlumistInds));

[totalBreadthPostP] = ranksum(totalBreadthPost(FluzoneInds),totalBreadthPost(FlumistInds));
[H3N2breadthPostP]  = ranksum(H3N2breadthPost(FluzoneInds),H3N2breadthPost(FlumistInds));


if(plotFlag)
  figInd = 1;

  [figInd] = plotSummaryResponsesByGroup(ptidRespMatPre(:,VicInds),groupNames([ptidStruct.preInd]),unique(groupNames([ptidStruct.preInd])),[projectName,' Pre-vaccination'],summaryMethod,figInd,yLimsSummary,fontSize);
  [figInd] = plotSummaryResponsesByGroup(ptidRespMatPost(:,VicInds),groupNames([ptidStruct.postInd]),unique(groupNames([ptidStruct.postInd])),[projectName,' Post-vaccination'],summaryMethod,figInd,yLimsSummary,fontSize);
  [figInd] = plotSummaryResponsesByGroup(ptidRespMatLT(:,VicInds),groupNames([ptidStruct.ltInd]),unique(groupNames([ptidStruct.ltInd])),[projectName,' LongTerm'],summaryMethod,figInd,yLimsSummary,fontSize);

  [figInd] = plotAllResponsesByGroup(ptidRespMatPre(:,VicInds),groupNames([ptidStruct.preInd]),unique(groupNames([ptidStruct.preInd])),[projectName,' Pre-vaccination'],figInd,yLims,fontSize);
  [figInd] = plotAllResponsesByGroup(ptidRespMatPost(:,VicInds),groupNames([ptidStruct.postInd]),unique(groupNames([ptidStruct.postInd])),[projectName,' Post-vaccination'],figInd,yLims,fontSize);
  [figInd] = plotAllResponsesByGroup(ptidRespMatLT(:,VicInds),groupNames([ptidStruct.ltInd]),unique(groupNames([ptidStruct.ltInd])),[projectName,' LongTerm'],figInd,yLims,fontSize);



  % plot total magnitude by treatment group
  figure(figInd);
  a = gca;
  set(a,'FontSize',16)
  boxplot(totalMagnitudePre,labels,'labels',treatmentLabels)
  title(sprintf('total magnitude Pre ranksum p = %.3f',totalMagPreP));
  figInd = figInd + 1; 

  figure(figInd);
  a = gca;
  set(a,'FontSize',16)
  boxplot(H3N2magnitudePre,labels,'labels',treatmentLabels);
  title(sprintf('H3N2 magnitude Pre ranksum p = %.3f',H3N2magPreP));
  figInd = figInd + 1; 

  % plot total magnitude by treatment group - post
  figure(figInd);
  a = gca;
  set(a,'FontSize',16)
  boxplot(totalMagnitudePost,labels,'labels',treatmentLabels)
  title(sprintf('total magnitude Post ranksum p = %.3f',totalMagPostP));
  figInd = figInd + 1; 

  figure(figInd)
  a = gca;
  set(a,'FontSize',16)
  boxplot(H3N2magnitudePost,labels,'labels',treatmentLabels);
  title(sprintf('H3N2 magnitude Post ranksum p = %.3f',H3N2magPostP));
  figInd = figInd + 1; 


  % plot total breadth by treatment group
  figure(figInd)
  a = gca;
  set(a,'FontSize',16)
  boxplot(totalBreadthPre,labels,'labels',treatmentLabels)
  title(sprintf('total breadth Pre ranksum p = %.3f',totalBreadthPreP));
  figInd = figInd + 1; 

  figure(figInd)
  a = gca;
  set(a,'FontSize',16)
  boxplot(H3N2breadthPre,labels,'labels',treatmentLabels)
  title(sprintf('total breadth Pre ranksum p = %.3f',H3N2breadthPreP));
  figInd = figInd + 1; 

  % plot total breadth by treatment group
  figure(figInd)
  a = gca;
  set(a,'FontSize',16)
  boxplot(totalBreadthPost,labels,'labels',treatmentLabels)
  title(sprintf('total breadth Post ranksum p = %.3f',totalBreadthPostP));
  figInd = figInd + 1; 

  figure(figInd)
  a = gca;
  set(a,'FontSize',16)
  boxplot(H3N2breadthPost,labels,'labels',treatmentLabels)
  title(sprintf('total breadth Post ranksum p = %.3f',H3N2breadthPostP));
  figInd = figInd + 1; 

  
  % PLOT dendrogram and corr matrices:
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStrPre  = 'Victoria H3N2 - Pre'; %dendrogram title
  titleStrPost = 'Victoria H3N2 - post'; %dendrogram title

  [figInd] = plotResponseDendrogram(ZstructVicPre,ptidLabels,figInd,fontSize,bwFlag,titleStrPre);
  [figInd] = plotResponseDendrogram(ZstructVicPost,ptidLabels,figInd,fontSize,bwFlag,titleStrPost);
  

  [Y1 I1] = sort(clusterLabelsVicPre);
  [figInd] = plotSampleCorrmat(figInd,corrMatVicPre(I1,I1),ptidLabels(I1),[-1 1],fontSize);
  title(titleStrPre);
  
  [Y2 I2] = sort(clusterLabelsVicPost);
  [figInd] = plotSampleCorrmat(figInd,corrMatVicPost(I2,I2),ptidLabels(I2),[-1 1],fontSize);
  title(titleStrPost);

  % plot correlation matrices of pre- and post- responses ordered by clusters induce by Pre-responses
  [figInd] = plotSampleCorrmat(figInd,corrMatVicPre(I1,I1),{ptidStruct(I1).ptid},[-1 1],fontSize);
  title(titleStrPre);

  [figInd] = plotSampleCorrmat(figInd,corrMatVicPost(I1,I1),{ptidStruct(I1).ptid},[-1 1],fontSize);
  title([titleStrPost, 'ordered by Pre']);
  
  [figInd] = plotAllResponsesByClusters(ptidRespMatPre(:,VicInds),clusterLabelsVicPre,'H3N2 Pre',figInd,yLims,fontSize)
  [figInd] = plotAllResponsesByClusters(ptidRespMatPost(:,VicInds),clusterLabelsVicPost,'H3N2 Post',figInd,yLims,fontSize)
  [figInd] = plotSummaryResponsesByClusters(ptidRespMatPost(:,VicInds),clusterLabelsVicPost,'H3N2 Post','median',figInd,yLimsSummary,fontSize);


  % when plotting POST clusters ordered by pre - plot the RAW post responses - not pre subtracted.
  ptidRawRespMat = reshape([ptidStruct.postResponse],size(ptidRespMatPost,2),[])';
  
  [figInd] = plotSummaryResponsesByClusters(ptidRawRespMat(:,VicInds),clusterLabelsVicPre,'H3N2 Post - ordered by Pre','median',figInd,yLimsSummary,fontSize);

  % for ERC application
  inds = find(clusterLabelsVicPre > 2); % first two clusters are very small.
  [figInd] = plotAllResponsesByClusters(ptidRespMatPre(inds,VicInds),clusterLabelsVicPre(inds),'',figInd,yLims,fontSize)

end


  %-------------------------------------------------------------------------------%  
  % plotting functions
  %-------------------------------------------------------------------------------%


   %print figures to files:
if (printFlag)
    figPath = [savePath,'FIG',filesep];
    for j=1:figInd
        figName = [figPath,projectName,'Figure_',num2str(j)];
        [figHandle] = printFigure(j,figName,{'eps','png'});
        close(figure(j));
    end
end



