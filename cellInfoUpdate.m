function Ps=cellInfoUpdate(Ps,I)

[xdim,ydim]=size(I);

for ci=1:1:numel(Ps)
    SingleContour=false(xdim,ydim);
    P=Ps{ci}.pts;    
    [R1,R2,NV]=getRibbon(P,Ps{ci}.thickness,xdim,ydim);
    contourList=[R1(:,:);R2(end:-1:1,:);R1(1,:)];
    for i=2:1:size(contourList,1)
        [xp,yp]=bresenham(contourList(i,1),contourList(i,2),contourList(i-1,1),contourList(i-1,2));
        pidx=sub2ind([xdim,ydim],xp,yp);
        SingleContour(pidx)=true;
    end
    SingleContour=imfill(SingleContour,'holes');
    
    Ps{ci}.intensity= mean(I(SingleContour>0));
    Ps{ci}.region = SingleContour;
    Ps{ci}.strip1=R1;
    Ps{ci}.strip2=R2;
    Ps{ci}.normvec = NV;
    Ps{ci}.length=sum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)));
end