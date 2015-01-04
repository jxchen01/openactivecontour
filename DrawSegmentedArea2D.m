function J=DrawSegmentedArea2D(Ps,I)
% Draw the contour
% (1) as regions in a logical image, 
% (2) as curves in the raw image (only display no return)
%
%  J=DrawSegmentedArea2D(P,I)
%
% inputs,
%  Ps : The conotours, N cells of size Ni x 2
%  I : The raw image
%
% outputs,
%  J : The binary image with all regions enclosed by contours
%
% Function is written by Jianxu Chen (University of Notre Dame) on Jan 2015
[xdim,ydim]=size(I);
J=false(xdim,ydim);

for i=1:1:numel(Ps)
    J(Ps{i}.region)=true;
    P=Ps{i}.pts;
    se=strel('disk',double(max([floor(Ps{i}.thickness),1])),0);
    tmp=false(xdim,ydim);
    tmp(round(P(1,1)),round(P(1,2)))=true;
    tmp(round(P(end,1)),round(P(end,2)))=true;
    J = J | imdilate(tmp,se);
end

figure(3);
imagesc(I,[0, 1]); axis off; axis equal; colormap(gray); 
hold on; contour(J,[1 1],'r','Linewidth',1); hold off;
drawnow



