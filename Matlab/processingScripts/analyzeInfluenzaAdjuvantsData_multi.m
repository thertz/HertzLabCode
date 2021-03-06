 
% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

clear all;
close all;

ZachFlag  = 1;
printFlag = 0;
plotFlag  = 0;
statFlag  = 0;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'Adjuvants';

if(ZachFlag)
 rootPath = [filesep,'\fhdata.fhcrc.org',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
else % Tomer path's
  rootPath = [filesep,'Volumes',filesep,'viddshared',filesep,'HertzLab',filesep]; % main rootPath under which data is located 
  savePath = [rootPath,'ArrayData',filesep,'Influenza',filesep,projectName,filesep];
  loadPath = savePath;
end
experimentDates = {'10_16_2014'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'rapa_zm_10_16_14_'}; % All slides and slideToSampleMapping file for each date have the exact same filePrefix
%09_12_2013 is rerunning of samples to make sure data is consistent
%experimentDates = {'07_23_2013','08_26_2013','08_30_2013','09_12_2013','09_03_2013'}; %dates of experiments, each one is a directory where GPR files are saved.
%expPrefix    = {'rapa_zm_07_23_13_','rapa_zm_08_26_13_','rapa_zm_08_30_13_','rapa_zm_09_12_13_','rapa_zm_09_03_13_'}; % All slides and slideToSampleMapping file for each date have the exact same filePrefix

 
fileSuffix = ''; % deprcated option to have a suffix to the name - we no longer allow this in the current naming format.
 
%-------------------------------------------------------------------------------%
% Additional parameters - project specific
%-------------------------------------------------------------------------------%
numArrays                = 2; %number of arrays on each slide.
typeFlag                 = 'median'; % uses the median over replicates of an antigen. Can also be 'mean'
colorFlags               = {'635'};
colorTags                = {'IgG'}; %which dye was used for which Ab type.
switchIDs_and_Names_Flag = 0; % old flag for problem with Scrips GPR files, deprecated
commentFlag              = 1; % slide to sample mapping contain a comment column. Current default, used for backward compatbility
fontSize                 = 16;
expGroups                = {'Sigma_Bc','Addavax_Bc','Alum_Bc','PBS_Bc','NC_Bc','Sigma_B6','Addavax_B6','Alum_B6','PBS_B6','NC_B6'};
expGroupsBcInds          = [1:5];  % which groups are BalbC in expGroups above! - HARD CODED
expGroupsB6Inds          = [6:10]; % which groups are B6 in expGroups above! - HARD CODED
adjuvantNames            = {'Sigma','Addavax','Alum','PBS','NC'};
yLims                    = [0 40000];
bgThreshold              = 5000; % threshold by which to flag antigens as ones that have high background and should
                                 % be removed from analysis (typically using a BSA negative control).
	
%samples to be included (1:14, 16:26, 29:60, 69:75, 77:81, 83, 85:94) - marked by Zach to remove duplicates
ptidInds = [1:14, 16:26, 29:60, 69:75, 77:81, 83, 85:94];


% Unsupervised clustering parameters:
numClusters          = 10; % number of clusters in the data (or expected number)
minResponseThreshold = 1000; %minimal threshold for responces to be included in clusterSamplesByResponseVectors, point below data is noise

%-------------------------------------------------------------------------------%

% HA peptide info
[pepData] = readHA_RapaPeptideData(rootPath);
%peptides in idices 1:96 in vectors have the following tube labels (important for matching!)
vietInds = strmatch('Influenza A/VietNam/1203/2004|H5N1|HA ',{pepData.strain});
x31Inds  = strmatch('Influenza A/aichi/2/68 (X31 recomb)|H3N2|HA',{pepData.strain});
PR8Inds  = strmatch('Influenza A PR8/ThomasLab|HA|H1N1',{pepData.strain});
%-------------------------------------------------------------------------------%
%defining elements outside of loop to add loop data too
 
slideNames      ={};
ptids           ={};
groupNames      ={};
mPeptResponses  =[];



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
end

% if this is set, limits the analysis to the specific ptids on this list, which can be from multiple experiment runs
if(~isempty(ptidInds))
  slideNames     = slideNames(ptidInds);
  ptids          = ptids(ptidInds);
  groupNames     = groupNames(ptidInds);
  mPeptResponses = mPeptResponses(ptidInds,:);
end

%Remove BSA background - accounting for any 2nd binding that is not related to the actual sample set
% note that this function can work an on array of arrayData structs, but here we are using it on a single arrayData struct.
[bgResponses] = computeBgResposnesFromArrayDataStruct(arrayData,'BSA'); % DONE OVER ENTIRE ARRAY USING ALL BSA samples
BSAinds = find(bgResponses.medianResp > bgThreshold);

%now zero out all responses in BSAinds across the entire dataset - high bg indices are considered non-reliable
% (conservative approach taken at this point)
mPeptResponsesBgSubtracted = mPeptResponses;
mPeptResponsesBgSubtracted(:,BSAinds) = 0;

%index all types of groups within
BcInds = [strmatch('Balb_c',ptids)]; 
B6Inds = [strmatch('Bk6',ptids) ];

Bc_NCinds = strmatch('NC_Bc',groupNames);
B6_NCinds = strmatch('NC_B6',groupNames); 
  
% Another sort of bg subtraction - currently highlighted as we found out most background comes from BSA, can later be
% merged with Bg subtraction
%subtract backgroung from each strain separately:
%mPeptResponsesBgSubtracted = mPeptResponses;
%mPeptResponsesBgSubtracted(BcInds,:) = mPeptResponses(BcInds,:) -repmat(median(mPeptResponses(Bc_NCinds,:)),length(BcInds),1);
%mPeptResponsesBgSubtracted(B6Inds,:) = mPeptResponses(B6Inds,:) -repmat(median(mPeptResponses(B6_NCinds,:)),length(B6Inds),1);


% cluster data using unserupvised correlation clustering and plot cluster dendrogram (tree) - need to define the
% following parameters up-top: numClusters and minThreshold
[distMat,Zstruct,clusterLabels] = clusterSamplesByResponseVectors(mPeptResponsesBgSubtracted,numClusters,minResponseThreshold);
[corrMat] = computeSampleCorrMatByResponseVectors(mPeptResponsesBgSubtracted,minResponseThreshold);




%statistically compare the responses of the two groups for each adjuvant separately:
if (statFlag)
  for i=1:length(adjuvantNames)
    
    currInds1 = strmatch([adjuvantNames{i},'_Bc'],groupNames); 
    currInds2 = strmatch([adjuvantNames{i},'_B6'],groupNames); 
    labels = [zeros(1,length(currInds1)) ones(1,length(currInds2))];
    
    groupResponses(i).name = adjuvantNames{i};  
    % ranksum test to compare the total responses across the array of two groups - in this case different strains of mice
    [groupResponses(i).pViet1203,~] = compareGroupsByTotalArrayResponses(mPeptResponsesBgSubtracted([currInds1;currInds2],vietInds),labels);
    
    %now filter responses blindly - first by percent responders (treatment blinded)
    [groupResponses(i).antigenFilters] = computeRespMatLabelBlindedFilters(mPeptResponsesBgSubtracted([currInds1;currInds2],vietInds),{'PercentResponders'},{0.2});
    [groupResponses(i).antigenSpecific_pValuesRankSum,groupResponses(i).antigenSpecific_pValuesFisher,groupResponses(i).group1responseHigher] = ...
	compareGroupsByAntigen(mPeptResponsesBgSubtracted([currInds1;currInds2],vietInds),labels,minResponseThreshold);

        
    %find peaks in maps:
    [groupResponses(i).peakResponsesBc] = findPeakResponses(mPeptResponsesBgSubtracted(currInds1,vietInds),5000);
    [groupResponses(i).peakResponsesB6] = findPeakResponses(mPeptResponsesBgSubtracted(currInds2,vietInds),5000);
  
    hotspotMap(i,:) = sum(groupResponses(i).peakResponsesBc) + sum(groupResponses(i).peakResponsesB6);
    
  end
  
  %now subtract out the negative control hotspots:
  ind = strmatch('NC',adjuvantNames);
  vaccInds = setdiff([1:length(adjuvantNames)],ind);
  hotspotMap(vaccInds,:) = hotspotMap(vaccInds,:) - repmat(hotspotMap(ind,:),length(vaccInds),1);
   
  %identify hotspots
  hotspotInds = find(sum(hotspotMap(vaccInds,:)) > 5);
  for i=1:length(hotspotInds)
    hotspotStruct(i).ind           = hotspotInds(i);
    hotspotStruct(i).seq           = pepData(hotspotInds(i)).sequence;
    hotspotStruct(i).begInd        = pepData(hotspotInds(i)).begInd;
  
  end
  
  fd = fopen([savePath,'sigAntigensAdjuvantStudy.txt'],'w')
  fprintf(fd,'Viet1203 peptides\n');
  fprintf(fd,'Sequence\tbegInd\t\n');
  for i=1:length(hotspotStruct)
    fprintf(fd,'%s\t%d\t\n',hotspotStruct(i).seq,hotspotStruct(i).begInd);
  end
  fclose(fd);
end






%----------------------------------------------------------------------------------------------------------%
%for plotting
if (plotFlag)
  figInd = 1;
  
  % plot all responses by group - limited only to the viet1203 inds
     
  % Plot all responses of each strain separately mice only
  [figInd] = plotAllResponsesByGroup(mPeptResponsesBgSubtracted(:,vietInds),groupNames, expGroups(expGroupsBcInds),'Adjuvants',figInd,yLims,fontSize); 
 
  % Plot all responses of each strain separately mice only
  [figInd] = plotAllResponsesByGroup(mPeptResponsesBgSubtracted(:,vietInds),groupNames, expGroups(expGroupsB6Inds),'Adjuvants',figInd,yLims,fontSize); 
  
  
  % compare responses for each group between the two strains
  for i=1:length(adjuvantNames)
    currInds = strmatch(adjuvantNames{i},expGroups);  
    [figInd] = plotAllResponsesByGroup(mPeptResponsesBgSubtracted,groupNames, expGroups(currInds),'Adjuvants',figInd,yLims,fontSize); 
  end
  
  %plot BSA respones (not BG subtracted)
  currInds = strmatch('BSA',groupNames);
  [figInd] = plotAllResponsesByGroup(mPeptResponses,groupNames,{'BSA','NC'},'Adjuvants',figInd,yLims,fontSize); 
  
  %[figInd] = plotAllResponsesByGroup(mPeptResponsesBgSubtracted(:,vietInds),groupNames,expGroups,'Adjuvants',figInd,yLims,fontSize);
  
  % PLOT dendrogram:
  bwFlag = 0; %print dendrograms in BW and not in blue.
  titleStr = ''; %dendrogram title
  fontSize = 8;
  [figInd] = plotResponseDendrogram(Zstruct,strrep(groupNames,'B6','Blk6'),figInd,fontSize,bwFlag,titleStr);
  
  %plot correlation matrix:    
  [Y I] = sort(clusterLabels);
  [figInd] = plotSampleCorrmat(figInd,corrMat(I,I),strrep(groupNames(I),'B6','Blk6'),[0.5 0.8],10);

  %plot hotstpots:
  for i=1:length(adjuvantNames)
    
    figure(figInd);
    subplot(5,1,i);
    bar(sum(groupResponses(i).peakResponsesBc));
    a = gca;
    set(a,'FontSize',16);
    set(a,'Ylim',[0 8]);
    title(['BALB/C ',groupResponses(i).name]);
    
    figure(figInd+1)
    subplot(5,1,i);
    bar(sum(groupResponses(i).peakResponsesB6));
    a = gca;
    set(a,'FontSize',16);
    set(a,'Ylim',[0 8]);
    title(['BLK6 ',groupResponses(i).name]);
  end
        
  
  %print figures to files:
  if (printFlag)
    figPath = [savePath,filesep,'FIG',filesep];
    for j=1:figInd
        figName = [figPath,projectName,'Figure_',num2str(j)];
        printFig(figure(j),figName,30,24);
    end
  end
end

