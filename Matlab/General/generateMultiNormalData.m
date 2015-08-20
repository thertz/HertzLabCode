% draws N points from a multi-variate gaussian of Ndim dimensions 
% with mean Mu (vector) and with Cov Matrix Sigma

function [data] = generateMultiNormalData(Mu,Sigma,N,Ndim)

%create the training set - first draw a random number of points
%then move them to the gaussians given using sigma and mu
data= randn(Ndim,N);

data= data + (Mu'*ones(1,N));
data= (data' * Sigma)';