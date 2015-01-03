function thickness=cellThickness(P,cellRegion,xdim,ydim)

ctl_idx=sub2ind([xdim,ydim],P(2:end-1,1),P(2:end-1,2));
tmp=zeros(xdim,ydim);
tmp(ctl_idx)=1;
distMap = bwdist(tmp);
    
bdList = bwboundaries(cellRegion);
bd=bdList{1};
bd_idx = sub2ind([xdim,ydim],bd(:,1),bd(:,2));
distValue = distMap(bd_idx);
thickness = mean(distValue(:));