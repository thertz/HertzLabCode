
% main script that loads up GRP file data and converts into arrayData struct which is a matlab representation of the
% GPR and experimental meta-data (ptids, slideNames etc.). This script should be used to create copies, each of which
% is for a given project. The script can load data from multiple experiments, so no need to create a new copy of this
% for every run, just for every project.

clear all;
close all;

%-------------------------------------------------------------------------------------------------------------------------------------%
% Parameters that must be set - Paths, dates etc. These need to be setup. As more dates are added simply add these to the list here
%-------------------------------------------------------------------------------------------------------------------------------------%
projectName = 'Rapa';


rootPath = [filesep,'Volumes',filesep,'data',filesep]; % main rootPath under which data is located 
savePath = [rootPath,'Work',filesep,'Results',filesep,'AntigenMicroarrays',filesep,projectName,filesep];
loadPath = [rootPath,'Work',filesep,'data',filesep,'AntigenMicroarrays',filesep,projectName,filesep];

experimentDates = {'11_06_2013','14_06_2013','25_06_2013','02_07_2013'}; %dates of experiments, each one is a directory where GPR files are saved.
expPrefix    = {'RAPA4_ZW_061113_','RAPA4_ZM_06_14_13_','Rapa_062513_ZW_','Rapa_ZW_070213_'}; % All slides and slideToSampleMapping file for each date have the exact same filePrefix

fileSuffix = ''; % deprcated option to have a suffix to the name - we no longer allow this in the current naming format.

%-------------------------------------------------------------------------------%
% Additional parameters - project specific
%-------------------------------------------------------------------------------%
numArrays                = 2; %number of arrays on each slide.
typeFlag                 = 'median'; % uses the median over replicates of an antigen. Can also be 'mean'
colorFlags               = {'635','532'};
colorTags                = {'IgG','IgM'}; %which dye was used for which Ab type.
switchIDs_and_Names_Flag = 0; % old flag for problem with Scrips GPR files, deprecated
commentFlag              = 1; % slide to sample mapping contain a comment column. Current default, used for backward compatbility
fontSize = 16;
%monoclonal contact sites:
[monoclonalStruct] = loadH5monoclonalContactSites();

% HA peptide info
[origArrayPepData] = readHA_RapaPeptideData('/users/thertz/');
%peptides in idices 1:96 in vectors have the following tube labels (important for matching!)
vietInds = strmatch('Influenza A/VietNam/1203/2004|H5N1|HA ',{origArrayPepData.strain});
x31Inds  = strmatch('Influenza A/aichi/2/68 (X31 recomb)|H3N2|HA',{origArrayPepData.strain});
PR8Inds  = strmatch('Influenza A PR8/ThomasLab|HA|H1N1',{origArrayPepData.strain});

monoNames = {monoclonalStruct.name};

for i=1:length(monoclonalStruct)
  for j=1:length(monoclonalStruct(i).contactSites)
    currContact = monoclonalStruct(i).contactSites(j).viet1203number;
    bInds = find(currContact >  [origArrayPepData(vietInds).begInd]);
    eInds = find(currContact <= [origArrayPepData(vietInds).endInd]);
    monoclonalStruct(i).arrayInds{j} = intersect(bInds,eInds);
  end
end


%-------------------------------------------------------------------------------%

for i=4:4%length(experimentDates)
  mappingFilename = [loadPath,experimentDates{i},filesep,expPrefix{i},'SlideToSampleMappings.txt'];
  currLoadPath    = [loadPath,experimentDates{i},filesep,'GPR',filesep];
  filePrefix      = [expPrefix{i},'slide_'];

  % name of file in which converted data will be stored for future use
  arrayDataFilename = [savePath,filesep,expPrefix{i},'arrayData_',typeFlag,'.mat']; 

  load(arrayDataFilename);

  slideNames = [arrayData.slideNames];
  ptids      = [arrayData.ptids];
  groupNames = [arrayData.groupNames];
  mPeptResponses      = [arrayData(1).responseMatrix{1}];

  if(i==3)
    baselineInds = [strmatch('NC',ptids);strmatch('NMS_20',ptids)];
  else
    baselineInds = [strmatch('NC',ptids);strmatch('NMS',ptids)];
  end
  
  for j=1:length(monoNames)
    allContactInds = unique([monoclonalStruct(j).arrayInds{:}]);
    
    if(i==3)
      currName = [monoNames{j}(1:4),'_20_',monoNames{j}(6:end)];
    else
      currName = monoNames{j};
    end
    monoInds = strmatch(currName,ptids);
    if(isempty(monoInds)) %VN04_2 instead of VN04-2
      monoInds = strmatch('VN04_2',ptids);
    end
    
    figure()
    bar(vietInds(allContactInds),mPeptResponses([monoInds;baselineInds],vietInds(allContactInds))');
    a = gca;
    set(a,'FontSize',fontSize);
    set(a,'yLim',[0 10000]);
    legend(strrep(ptids([monoInds;baselineInds]),'_',' '));
    
    figName = [savePath,experimentDates{i},filesep,expPrefix{i},monoNames{j},'_Responses'];
    printFig(gcf,figName,20,16);
  end
  
end

