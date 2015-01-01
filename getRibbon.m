function [R1,R2]=getRibbon(P,thickness,xdim,ydim)

N=GetContourNormals2D(P); %%% unit normal vectors
R1=N.*thickness;
R2=N.*(-thickness);

R1(:,1)=min(max(R1(:,1),1), xdim);
R1(:,2)=min(max(R1(:,2),1), ydim);

R2(:,1)=min(max(R2(:,1),1), xdim);
R2(:,2)=min(max(R2(:,2),1), ydim);
