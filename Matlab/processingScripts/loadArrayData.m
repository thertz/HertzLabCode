
function [arrayData,sampleInds] = loadArrayData(mappingFileName,filePrefix,loadPath,processGPRfunctionHandle,typeFlag,colorFlag,switchIDs_and_Names_Flag)

%load array data from a given list of slides. assumes by convention that the slideNames always end with slide number
%followed by a .gpr extension.
% params:
% mappingFileName          - name of file in conventional format that maps samples to slides and epxerimental conditions(groups)
% filePrefix               - prefix of array file Names
% loadPath                 - path to GPR files, including rootPath
% processGPRfunctionHandle - function handle to processing of the raw GPR file, which is different for each project.
% typeFlag                 - can be mean/median to indicate how to handle replicates of each antigen.
% colorFlag                - which channel to work on, 532 or 635
% switchIDs_and_Names_Flag - some GPR files based on scripps printing have the IDs and Names fields switched. If
% raised, will swap these for all code to work properly.

% returns:
% sampleInds - list of gpr files found and used for analysis.

if(~exist('switchIDs_and_Names_Flag','var'))
  switchIDs_and_Names_Flag = 0;
end

[slideNames,ptids,groupNames] = parseSlideMappingFile(mappingFileName);

currInd = 1; %used ot index existing GPR files.
remInds = [];
for s=1:length(slideNames)
  
    gprFilename = [loadPath,filePrefix,slideNames{s},fileSuffix{s},'.gpr'];

  try
    [gprStruct] = gprread(gprFilename,'cleanColNames',1);

    if(switchIDs_and_Names_Flag)
      [gprStruct] = switchGPR_IDs_and_Names(gprStruct);
      [gprStruct] = cleanGPRstructScrippsVersion(gprStruct);
    else
      [gprStruct] = cleanGPRstruct(gprStruct);
    end
    
    [seqGPRstruct,otherGPRstruct] = processGPRfunctionHandle(gprStruct);
    
    [meanResponses,antigenNames] = computeMeanResponsesByNames(seqGPRstruct,typeFlag,colorFlag);
    
    arrayData.seqGPRstructs{currInd}   = seqGPRstruct;
    arrayData.otherGPRstructs{currInd} = otherGPRstruct;
    
    if(currInd==1)
      arrayData.antigenNames = antigenNames;
    else %check to verify that antigenNames are consistent across all slides!
      if(sum(strcmp(antigenNames,arrayData.antigenNames)) ~= length(arrayData.antigenNames))
	error(sprintf('antigen names of slide %s are inconsistent, exiting',gprFilename));
     end
    end
    
    arrayData.responseMatrix(currInd,:) = meanResponses;
    currInd = currInd + 1;
  catch
    
    disp(sprintf('gpr File: %s not found',gprFilename));
    remInds = [remInds s];
  end
end

sampleInds = setdiff([1:length(slideNames)],remInds);

slideNames(remInds) = [];
ptids(remInds)      = [];
groupNames(remInds) = [];

arrayData.slideNames = slideNames';
arrayData.ptids      = ptids';
arrayData.groupNames = groupNames';
