%get shape context descriptor histograms (Belongie, Malik)
%- adapted from michael fink's code. 
% this version samples the descriptor using an n x n  defined by descX
% descY grid.

% input:
% edgeX,edgeY        - indices of all edge pixles in the image
% descX,descY        - indices of descriptors - on a uniform grid.
% numOrientationBins - number of orientation bins. default set to 12.
% numRadiusBins      - number of radius bins. defauls set to 5.

% outputs:
% scHistogram - the Shape Context Histograms. Each histogram is a column
% in the matrix. matrix of size (numOrientationBins*numRadiusBins) x (numFeatures)
function [ShapeContextHistograms]=getShapeContextDescriptors(edgeX,edgeY,descX,descY,numOrientationBins,numRadiusBins)

if(~exist('numOrientationBins','var'))
  numOrientationBins= 12;
end

if(~exist('numRadiusBins','var'))
  numRadiusBins= 5;
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

numFeatures   = length(descX);

ShapeContextHistograms = zeros(numOrientationBins*numRadiusBins,numFeatures); 

for i=1:numFeatures
    tx = edgeX-descX(i);
    ty = edgeY-descY(i);
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