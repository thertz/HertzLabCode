

%filter a dataset according to some filtering method - could be classSize,
%or supertype, or whatever. FilteringParams is a struct with various
%filtering parameters required for various filtering methods...
function [peptides,labels,labelMap] = filterDataset(origPeptides,origLabels,origLabelMap,filteringMethod,FilteringParams)

allInds = [];
bigInds =[];
switch(filteringMethod)
    
    case 'classSize'
        
        %select only based on the number of positive examples of course!
        plg     = pilug(origLabels);
        posInds = find(plg(:,1) > 0);
        plg     = plg(posInds,:);
        
        bigInds = find(plg(:,2) >= FilteringParams.classSize);
        for i=1:length(bigInds)
            posInds  = find(origLabels == plg(bigInds(i),1));
            negInds  = find(origLabels == -plg(bigInds(i),1));
            allInds = [allInds ;posInds ;negInds];
        end
    case 'supertype'
       %A2.alleles = {'A0201','A6802','A0202','A6901','A0203','A0204','A0205','A0206','A0207','A0209'};
       A2.alleles = {'A02','A6802','A6901'};
        for i=1:length(A2.alleles)
            currInd = strmatch(A2.alleles{i},origLabelMap,'exact');
            if(isempty(currInd))
                currInd = strmatch(A2.alleles{i}(1:3),origLabelMap,'exact');
            end
            
            bigInds = [bigInds ; currInd];
            posInds  = find(origLabels == currInd);
            negInds  = find(origLabels == -currInd);
            allInds = [allInds ;posInds ;negInds];
        end
        
    case 'singleAllele'
        currInd = strmatch(FilteringParams.allele,origLabelMap,'exact');
        posInds  = find(origLabels == currInd);
        negInds  = find(origLabels == -currInd);
        allInds = [allInds ;posInds ;negInds];
    
    case 'alleleList'
        for i=1:length(FilteringParams.alleleList)
            currInd = strmatch(FilteringParams.alleleList{i},origLabelMap,'exact');
            bigInds = [bigInds; currInd];
            posInds  = find(origLabels == currInd);
            negInds  = find(origLabels == -currInd);
            allInds = [allInds ;posInds ;negInds];
        end
end

peptides = origPeptides(allInds,:);
labels   = origLabels(allInds);
labelMap = origLabelMap(bigInds);
