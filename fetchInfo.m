function [RepelForce, ContourImage]=fetchInfo(Ps,BMap)

%%% key parameter %%%
repelThresh=5;
%%%%%%%%%%%%%%%%%%%%%

[xdim,ydim]=size(BMap);
numContour=numel(Ps);

RepelForce = zeros(xdim,ydim,numContour+1);
RepelForce(:,:,end)=BMap(:,:);

% build the augmented regions
ContourImage=false(xdim,ydim);
for ci=1:1:numContour
    ContourImage = ContourImage | Ps{ci}.region ;
    
    repelMap = bwdist(Ps{ci}.region);
    repelMap(repelMap>repelThresh)=0;
    idx=find(repelMap>0);
    repelMap(idx)=1./(1+exp(3.*(repelMap(idx)-3)));
    repelMap(Ps{ci}.region)=100;
    RepelForce(:,:,ci)=repelMap(:,:);
end