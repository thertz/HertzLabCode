%clutsers the data using the required method. is a wrapper for Linkage from matlab

%params:
% M - matrix that can be a similarityMatrix or a distance matrix.
% Mtype - the matrix type - 'similarity' or 'distance'
% k - number of required clusters.
% method - clustering method: complete - complete linkage
%                             single   - single linkage.
%                             average - average linkage
%                             ward - ward's method
%                             all - all methods above!
% labels - the data labels (optional)

% returns:
% parition     - a vector of data size which contains the cluster assignment of each data point.
% clusterStats - purity accuracy and Z-score in a 3xmethods vector.
% Z            - linkage format in matlab that allows to produce dendrograms.
function [clusterStats,partition,Z]= LinkageClustering(M,Mtype,k,method,labels)

  if(~exist('labels','var'))
    labels= [];
  end
  
  if(strcmp(Mtype,'similarity'))
    %turn simlrity into distance:
    %begin by moving all similarities to pos numbers: not for all sim to dist transformations.
    tempMat= M+ abs(min(min(M))) + 1;
    tempMat= tempMat./max(max(tempMat));
    
    % %few options here:
    distMat= 1 - tempMat;
    %distMat= 1./tempMat;
    %distMat= (max(max(tempMat)) - tempMat) + 1;
    %distMat= log(1./tempMat);
    
    N= size(distMat,1);
    distMat([1:N+1:N^2])=0;
    
  else
    distMat= M;
  end
  
  if(strcmp(method,'all'))
    methods= {'single','average','complete','ward'};
  else
    methods= {method};
  end
  
  for i=1:length(methods)
    [clusterStats(i,:), partition(i,:), Z{i}]= matCluster(distMat,k, methods{i},labels);
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clusterStats,clusters,Z]= matCluster(distMat, k, method, labels)
  
%transfom into pdist format:
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

