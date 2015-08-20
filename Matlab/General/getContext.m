%get shape context descriptor histograms (Belongie, Malik)
%- adapted from michael fink's code.

% input:
% x,y                - indices of all edge pixles in the image
% numFeatures        - number of descriptors required. usually between 30-100.
% numOrientationBins - number of orientation bins. default set to 12.
% numRadiusBins      - number of radius bins. defauls set to 5.

% outputs:
% scHistogram - the Shape Context Histograms. Each histogram is a column
% in the matrix. matrix of size (numOrientationBins*numRadiusBins) x (numFeatures)
function [ShapeConstextHistograms,edgeInds]=getShapeContextDescriptors(edgeX,edgeY,numFeatures,numOrientationBins,numRadiusBins)

if(~exist('numOrientationBins','var'))
  numOrientationBins= 12;
end

if(~exist('numRadiusBins','var'))
  numOrientationBins= 5;
end

% USER SELECTS POINTS????
% if nEdge<0 % user select points    
%     nEdge = -nEdge;
%     nInd = zeros(1,nEdge);
%     for iEdge=1:nEdge
%         [ux,uy,b] = ginput(1);
%         [nnDist nnInd] = min((x-ux).^2 + (y-uy).^2);
%         nInd(iEdge) = nnInd;
%         for iT=1:nT
%             TT = iT*2*pi/nT;
%             line(x(nnInd)+[0 cos(TT)]*(2^nR-1),y(nnInd)+[0 sin(TT)]*(2^nR-1),'color',[1 0 0],'lineWidth',2);
%         end
%         for iR=1:nR
%             RR = 2^iR-1;
%             rectangle('Position', [x(nnInd)-RR y(nnInd)-RR 2*RR 2*RR],'Curvature',[1 1],'EdgeColor',[1 0 0],'lineWidth',2);
%         end
%     end
%else

%choose edge pixels on which to build descriptors:

numEdgePixels = length(edgeX); % number of edge pixels in the image.
edgeInds      = randperm(numEdgePixels);
numFeatures   = min(numFeatures,numEdgePixels);
edgeInds = edgeInds(1:numFeatures); 

ShapeContextHistograms = zeros(numOrientationBins*numRadiusBins,numFeatures); 

for i=1:numFeatures
    tx = x-x(edgeInds(i));
    ty = y-y(edgeInds(i));
    tx = [tx(1:edgeInds(i)-1);tx(edgeINds(i)+1:end)];
    ty = [ty(1:edgeInds(i)-1);ty(edgeInds(i)+1:end)];
    tR = ceil(log2(1+sqrt(tx.^2+ty.^2)));
    tT = ceil(numOrientationBins/2*(atan2(tx,ty)/pi+1));
    tA = zeros(numOrientationBins,numRadiusBins);

    for iT=1:numOrientationBins
        for iR=1:numRadiusBins
            tA(iT,iR) = sum(tR==iR & tT==iT)/(2^(iR/2));
        end
    end
    tA = tA-mean(tA(:));    tA = tA-std(tA(:)); 
    ShapeContextHistograms(:,i)=reshape(tA,[],1);
end