%compute shape context descriptors using a fixed grid of points inputed
%on the input image.
function [SChist] = computeShapeContext(Image,GridPoints);
  

%adapted from Malik code:
  
  nGrid = size(GridPoints,1);
  
  mean_dist_global = [];
  nbins_theta      = 12;
  nbins_r          = 5;
  r_inner          = 1/8;
  r_outer          = 2;

  %disp('extracting boundary points...')
  [x1,y1,t1]=bdry_extract_3(Image);  

  BoundaryPoints = [x1,y1];
  nBoundary      = size(BoundaryPoints,1);
  out_vec        = zeros(1,nBoundary);
    
  [SChist,meanDist] = sc_compute1(GridPoints', BoundaryPoints' ,mean_dist_global,nbins_theta,nbins_r,r_inner,r_outer,out_vec);
  

