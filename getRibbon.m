function [R1,R2,N]=getRibbon(P,thickness,xdim,ydim)

N=GetContourNormals2D(P); %%% unit normal vectors
dt=N.*thickness;
R1=P+dt;
R2=P-dt;

R1(:,1)=min(max(R1(:,1),1), xdim);
R1(:,2)=min(max(R1(:,2),1), ydim);

R2(:,1)=min(max(R2(:,1),1), xdim);
R2(:,2)=min(max(R2(:,2),1), ydim);
