function N=GetContourNormals2D(P)
% This function calculates the normals, of the contour points
% using the neighbouring points of each contour point
%
% N=GetContourNormals2D(P)
% 
% inputs,
%  P : List with contour coordinates M x 2
%
% outputs,
%  N : List with contour normals M x 2
%
% Function is written by D.Kroon University of Twente (July 2010)
% Modified by Jianxu Chen (University of Notre Dame) at Jan 2015

% Use the a'th neighbour to calculate the normal (more stable)
a=3;

% From array to separate x,y
xt=P(:,1); yt=P(:,2);

% Derivatives of contour
n=length(xt);
dx=zeros(n,1);
dy=zeros(n,1);

hh=2*a;
dx(a+1:end-a)=(xt(2*a+1:end)-xt(1:end-2*a))./hh;
dy(a+1:end-a)=(yt(2*a+1:end)-yt(1:end-2*a))./hh;

for i=a:-1:2
    hh=(2*i-2);
    dx(i)=(xt(i+i-1)-xt(1))/hh; dx(end-i+1)=(xt(end)-xt(end-2*i+2))/hh;
    dy(i)=(yt(i+i-1)-yt(1))/hh; dy(end-i+1)=(yt(end)-yt(end-2*i+2))/hh;
end

dx(1)=xt(2)-xt(1);dx(end)=xt(end)-xt(end-1);
dy(1)=yt(2)-yt(1);dy(end)=yt(end)-yt(end-1);

% Normals of contourpoints
l=sqrt(dx.^2+dy.^2);
nx = -dy./l; 
ny =  dx./l;
N(:,1)=nx; N(:,2)=ny;
