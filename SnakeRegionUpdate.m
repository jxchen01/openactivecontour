function P=SnakeRegionUpdate(I,B,P,gamma,kappa,delta,thickness,targetLength)
% This function will calculate one iteration of ribbon Snake movement
%
% P=SnakeRegionUpdate(I,B,P,gamma,kappa,thickness,targetLength)
%
% inputs,
%   I : Raw Image in grayscale
%   B : Internal force (smoothness) matrix
%   P : The contour points N x 2;
%   gamma : Time step
%   kappa : External (image) field weight
%   delta : stretching force due to length prior
%
% outputs,
%   P : The new (moved) contour points N x 2;
%
% Function is written by Jianxu Chen (University of Notre dame) on Jan 2015
% modified based on a script from D.Kroon University of Twente 

[xdim,ydim]=size(I);
% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),xdim);
P(:,2)=min(max(P(:,2),1),ydim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate force by Chan-Vese on the augmented region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ContourImage=zeros(xdim,ydim);

% get the contour region
[R1,R2,NV]=getRibbon(P,thickness,xdim,ydim);
contourList=[R1(:,:);R2(end:-1:1,:);R1(1,:)];
for i=2:1:size(contourList,1)
    [xp,yp]=bresenham(contourList(i,1),contourList(i,2),contourList(i-1,1),contourList(i-1,2));
    pidx=sub2ind([xdim,ydim],xp,yp);
    ContourImage(pidx)=1;
end
clear i xp yp pidx

se1=strel('disk',1,0);
ContourImage = imclose(ContourImage,se1);
ContourImage=imfill(ContourImage);

interiorIntensity = mean(I(ContourImage>0));
exteriorIntensity = mean(I(ContourImage==0));

dphi_all = (I-interiorIntensity).^2-(I-exteriorIntensity).^2;
dphi(:,1) = interp2(dphi_all,R1(:,2),R1(:,1));
dphi(:,2) = interp2(dphi_all,R2(:,2),R1(:,1));
% Interp2, can give nan's if contour close to border
dphi(isnan(dphi))=0;

F_CV_Norm = -(dphi(:,1)-dphi(:,2))./max(abs(dphi(:)));
F_CV(:,1)=NV(:,1).*F_CV_Norm;
F_CV(:,2)=NV(:,2).*F_CV_Norm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the stretching force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len=sum(sqrt(sum((P(2:end,:)-P(1:end-1,:)).^2,2)));
if(len<0.85*targetLength)
    phi=1;
elseif(len<0.95*targetLength)
    phi=0.5;
elseif(len>1.05*targetLength)
    phi=-1;
else % expected length [0.95L, 1.05L]
    phi=0;
end
nn1=P(1,:)-P(2,:);
sf1=(delta*phi).*(nn1./hypot(nn1(1),nn1(2)));
nn2=P(end,:)-P(end-1,:);
sf2=(delta*phi).*(nn2./hypot(nn2(1),nn2(2)));

% Update contour positions
ssx = gamma*P(:,1) + F_CV(:,1);% + Fext2(:,1);
ssy = gamma*P(:,2) + F_CV(:,2);% + Fext2(:,2);
ssx(1) = ssx(1)+sf1(1); ssy(1) = ssy(1)+sf1(2);
ssx(end) = ssx(end)+sf2(1); ssy(end) = ssy(end)+sf2(2);
P(:,1) = B * ssx;
P(:,2) = B * ssy;

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),xdim);
P(:,2)=min(max(P(:,2),1),ydim);

    
