%computes PCA using a trick to overcome large matrices where the number
%of dimensions of each data point is larger than the number of data points
%the resulting eigenvectors are orthonormalized.

% params:   TrainingSet: data matrix, each data point is a column.
%           numPC: number of PC's to project on.

% returns : ReducedTrainingSet: projected data
%           meanI - the mean data vector (For recostruction)
%           W - the eigenvectors used for projection

function [ReducedTrainingSet,W,meanI]= PCA(TrainingSet,numPC)

[rows cols]= size(TrainingSet);

meanI= mean(TrainingSet');

%determine if to work on original feature space (Columns) or on data
%space (rows)...
transPCAflag = ( (rows > cols) && (rows > 2000) );

%if feature set is larger than dataset size...
if(transPCAflag)
  for i=1:cols
    i;
    TrainingSet(:,i)= TrainingSet(:,i) - meanI';
  end
else
  TrainingSet= TrainingSet - meanI'*ones(1,cols);
end

if(transPCAflag)
  temp      = TrainingSet'*TrainingSet;
  [U Sig V] = svd(temp);

  U = TrainingSet*U;  
else
  temp      = TrainingSet*TrainingSet';
  [U Sig V] = svd(temp);
end

%normalize the eigenvectors :
for i=1:numPC
    U(:,i)= U(:,i)./norm(U(:,i));
end

W = U(:,1:numPC);

%project the data on the eigenvectors:
ReducedTrainingSet = (W'*TrainingSet);





