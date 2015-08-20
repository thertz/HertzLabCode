
function [arrayData,sampleInds] = loadMultiplexArrayData(mappingFilename,filePrefix,fileSuffix,loadPath,numArrays,processGPRfunctionHandle,typeFlag,colorFlags,switchIDs_and_Names_Flag,commentFlag)

%load array data from a given list of slides printed with multiple arrays. assumes by convention that the slideNames always end with slide number
%followed by a .gpr extension.
% params:
% mappingFilename          - name of file in conventional format that maps samples to slides and epxerimental conditions(groups)
% filePrefix               - prefix of array file Names
% loadPath                 - path to GPR files, including rootPath
% processGPRfunctionHandle - function handle to processing of the raw GPR file, which is different for each project.
% typeFlag                 - can be mean/median to indicate how to handle replicates of each antigen.
% colorFlag                - which channel to work on, 532 or 635
% switchIDs_and_Names_Flag - some GPR files based on scripps printing have the IDs and Names fields switched. If
% raised, will swap these for all code to work properly.
% commentFlag - if up, mapping file has a comment column (new as of May 2013)

% returns:
% sampleInds - list of gpr files found and used for analysis.
% 
% Modified Dec. 17th 2014 to include block numbers into struct for local background normalization
if(~exist('switchIDs_and_Names_Flag','var'))
  switchIDs_and_Names_Flag = 0;
end

if(~exist('commentFlag','var'))
  commentFlag = 0;
end


[slideNames,ptids,groupNames,arrayNumbers,comments] = parseMultiplexSlideMappingFile(mappingFilename,commentFlag);

currInd = 1; %used ot index existing GPR files.
remInds = [];
for s=1:length(slideNames)
  
    gprFilename = [loadPath,filePrefix,slideNames{s},fileSuffix,'.gpr'];

  try
    [gprStruct] = gprread(gprFilename,'cleanColNames',1);

    if(switchIDs_and_Names_Flag)
      [gprStruct] = switchGPR_IDs_and_Names(gprStruct);
      [gprStruct] = cleanGPRstructScrippsVersion(gprStruct);
    else
      [gprStruct] = cleanGPRstruct(gprStruct);
    end

    [seqGPRstruct,otherGPRstruct] = processGPRfunctionHandle(gprStruct,arrayNumbers(s),numArrays);

    %THIS IS WHERE I'M AT - HANDLE CASE OF OTHER GPR STRUCT EMPTY - AND GO BACK TO FLU HA RAPA TO HANDLE THIS CASE TO
    %PREVENT CODE BREAKING HERE.    
    arrayData.seqGPRstructs{currInd}   = seqGPRstruct;
    arrayData.otherGPRstruct{currInd}  = otherGPRstruct;

    
    otherAntigenNames = [];
    if(length(colorFlags)==1)
      [meanResponses,antigenNames,blockNumbers] = computeMeanResponsesByNames(seqGPRstruct,typeFlag,colorFlags{1});
      arrayData.responseMatrix{1}(currInd,:) = meanResponses;
      
      if(~isempty(otherGPRstruct))
       [meanResponses,otherAntigenNames,otherBlockNumbers] = computeMeanResponsesByNames(otherGPRstruct,typeFlag,colorFlags{1});
       arrayData.otherResponseMatrix{1}(currInd,:) = meanResponses;
      end
    else
      for c=1:length(colorFlags)
       [meanResponses,antigenNames,blockNumbers] = computeMeanResponsesByNames(seqGPRstruct,typeFlag,colorFlags{c});
       arrayData.responseMatrix{c}(currInd,:) = meanResponses;

       if(~isempty(otherGPRstruct))
         [meanResponses,otherAntigenNames,otherBlockNumbers] = computeMeanResponsesByNames(otherGPRstruct,typeFlag,colorFlags{c});
         arrayData.otherResponseMatrix{c}(currInd,:) = meanResponses;
       end
      end
    end
   
    %load some labels only on first array
    if(currInd==1)
      arrayData.antigenNames      = antigenNames;
      arrayData.otherAntigenNames = otherAntigenNames;
      arrayData.colorFlags        = colorFlags;
      arrayData.blockNumbers      = blockNumbers;
      if(~isempty(otherGPRstruct))
        arrayData.otherBlockNumbers = otherBlockNumbers;
      end
    else %check to verify that antigenNames are consistent across all slides!
      if(sum(strcmp(antigenNames,arrayData.antigenNames)) ~= length(arrayData.antigenNames))
       error(sprintf('antigen names of slide %s are inconsistent, exiting',gprFilename));
     end
   end
    
   currInd = currInd + 1;
  catch
    disp(sprintf('gpr File: %s not found',gprFilename));
    keyboard;
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
arrayData.comments   = comments;

