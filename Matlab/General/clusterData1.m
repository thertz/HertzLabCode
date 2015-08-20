%clutsers the data using the required method. is a wrapper for Linkage from matlab

%params:
% M - matrix that can be a similarityMatrix or a distance matrix.
% Mtype - the matrix type - 'similarity' or 'distance'
% k - number of required clusters.
% method - clustering method: complete - complete linkage
%                             single   - single linkage.
%                             average - average linkage
%                             ward - ward's method
%                             spc - spc clustering.
%                             yoram - typical cut.
%                             all - all methods above!
% labels - the data labels.

% returns:
% parition - a vector of data size which contains the cluster assignment of each data point.
% clusterStats - purity accuracy and Z-score in a 3xmethods vector.
function [clusterStats,partition]= clusterData1(M,Mtype,k,labels,method,dataName)

if(strcmp(Mtype,'similarity'))
  %turn simlrity into distance:

  %begin by moving all similarities to pos numbers: not for all sim to dist transformations.
  tempMat= M+ abs(min(min(M))) + 1;

  tempMat= tempMat./max(max(tempMat));

  distMat= 1 - tempMat;
  
  %distMat= max(max(tempMat)) - tempMat;
  %distMat(find(distMat==0))=eps;
  
  % %few options here:
  %distMat= 1./tempMat;
  %distMat= (max(max(tempMat)) - tempMat) + 1;
  %distMat= log(1./tempMat);
  N= size(distMat,1);
  distMat([1:N+1:N^2])=0;

else
  distMat= M;
  %distMat= distMat./max(max(distMat));
end

if(strcmp(method,'all'))
  methods= {'single','average','complete','ward'};
else
  methods= {method};
end

for i=1:length(methods)
  
  if(strcmp(methods{i},'spc'))
    disp('run spc')
    Tmin=1;
    deltaT= 0.1;
    Tmax= 1;
    numNeighbors=8;
    dataName= ['data',num2str(length(labels))]

    if(strcmp(Mtype,'similarity'))
      [SPCclusterStats, SPCclusters]= SPCcluster(tempMat,numNeighbors,labels,dataName,Tmin,Tmax,deltaT);   
    end

  elseif (strcmp(methods{i},'yoram'))

    %tweak for each DATA sep: (voodoo!)
    if(strcmp(dataName,'proteinData'))
      numNeighbors= [3:5];
      mkMethod='threshold';
      
    elseif(strcmp(dataName,'ionosphereData'))
      numNeighbors=[3:5];
      mkMethod='threshold';
    elseif(strcmp(dataName,'balanceData'))
      numNeighbors=[8:10];
      mkMethod='knn';
    elseif(strcmp(dataName,'bostonData'))
      numNeighbors=[3:8];
      mkMethod='threshold';
    elseif(strcmp(dataName,'pimaDiabetesData'))
      numNeighbors=[3:5];
      mkMethod='threshold';
    elseif(strcmp(dataName,'breastData'))
      numNeighbors=[5];
      mkMethod='threshold';
    end
      %actually passes the similarity matrix.
    [clusterStats, partition]= TypicalCluster(M,numNeighbors,labels,dataName,mkMethod);

  else
    [clusterStats(i,:), partition(i,:)]= matCluster(distMat,k, labels, methods{i});
  end
end

%if only SPC use different temparatures instead of different methods - just for now!
if(strcmp(method,'spc'))
  clusterStats= SPCclusterStats;
  partition= SPCclusters;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clusterStats,clusters]= matCluster(distMat, k, labels, method)

%transfomr into pdist format:
%nice trick to overcome zeros in distances:
zeroInds= find(distMat==0);
distMat(zeroInds)=nan;

tempMat= triu(distMat,1);
tempMat= tempMat';
Y= tempMat(:)';
Y(find(Y==0))=[];
Y(isnan(Y))= 0;

Z = linkage(Y,method);

clusters = cluster(Z, k);
clusters= clusters';

if(~isempty(labels))
  [pu , ac]=purity_accuracy2(labels,clusters,k);
  if(pu+ac == 0)
    zscore=0;
  else
    zscore=2*pu*ac/(pu+ac);
  end
else
  pu=-1;
  ac=-1;
  zscore=-1;
    
end
clusterStats(1,1)=pu;
clusterStats(1,2)=ac;
clusterStats(1,3)=zscore;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

