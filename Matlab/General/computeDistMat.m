%computes mahalnobis distance for data:

function [distMat]= computeDistMat(data,A)

n= size(data,1);

% prepare distance matrix
distances=zeros(n,n);
for i=1:n
     
      temp=ones(n,1)*data(i,:);
      x_minus_y=temp-data;
      distances(i,:)= sum( ( ( x_minus_y*A ).*x_minus_y )' );
end

distMat= distances;