function P=SnakeRegionUpdate(I,B,Ps,BMap,gamma,kappa,delta)
% This function will calculate one iteration of ribbon Snake movement
%
% P=SnakeRegionUpdate(I,B,P,BMap,gamma,kappa)
%
% inputs,
%   I : Raw Image in grayscale
%   B : Internal force (smoothness) matrix
%   Ps : The contours, N cells with array 'pts' of size Ns x 2;
%   BMap : barrier maps from matched cells and other moving contours
%   gamma : Time step
%   kappa : weight of repelling force
%   delta : weight of stretching force due to length prior
%
% outputs,
%   Ps : The new (moved) contours
%
% Function is written by Jianxu Chen (University of Notre dame) on Jan 2015
% modified based on a script from D.Kroon University of Twente 

[xdim,ydim]=size(I);

% get repel force, cell region (filled region, two strips, normal vector)
% and interior intensity of each cell
numContour=numel(Ps);
nPoints=size(Ps{1}.pts,1);

[RepelForce,ContourImage] = fetchInfo(Ps,BMap);

% get the exterior intensity
exteriorIntensity = mean(I(~ContourImage));

% apply forces on each cell and update point list
dphi=zeros(nPoints,2);
rp=zeros(nPoints,2);
Fext=zeros(nPoints,2);
for ci=1:1:numContour
    % retrieve info of current cell
    interiorIntensity = Ps{ci}.intensity;
    P =Ps{ci}.pts;
    R1=Ps{ci}.strip1;
    R2=Ps{ci}.strip2;
    NV=Ps{ci}.normvec;
    len=Ps{ci}.length;
    targetLength = Ps{ci}.targetLength;
        
    % Calculate image force based on Chan-Vese model
    dphi_all = (I-interiorIntensity).^2-(I-exteriorIntensity).^2;
    dphi(:,1) = interp2(dphi_all,R1(:,2),R1(:,1));
    dphi(:,2) = interp2(dphi_all,R2(:,2),R1(:,1));    
    dphi(isnan(dphi))=0;
    
    F_CV_Norm = -(dphi(:,1)-dphi(:,2))./max(abs(dphi(:)));
    
    % Calculate repeling force
    id=setdiff(1:1:numContour+1,ci);
    rp_all=sum(RepelForce(:,:,id),3);
    rp(:,1) = interp2(rp_all,R1(:,2),R1(:,1));
    rp(:,2) = interp2(rp_all,R2(:,2),R1(:,1)); 
    rp(isnan(rp))=0;
    
    F_repel_Norm = -(rp(:,1)-rp(:,2));
    
    % Combing repeling force and image force into External Force Vector
    Fext(:,1)=NV(:,1).*(F_CV_Norm+kappa.*F_repel_Norm);
    Fext(:,2)=NV(:,2).*(F_CV_Norm+kappa.*F_repel_Norm);
    
    
    % Calculate the stretching force
    if(len<0.85*targetLength)
        phi=1;
    elseif(len<0.99*targetLength)
        phi=0.5;
    elseif(len>1.1*targetLength)
        phi=-1;
    else % expected length [0.95L, 1.05L]
        phi=0;
    end
    % stretch on head
    nn1=P(1,:)-P(2,:); sf1=(delta*phi).*(nn1./hypot(nn1(1),nn1(2)));
    % stretch on tail
    nn2=P(end,:)-P(end-1,:); sf2=(delta*phi).*(nn2./hypot(nn2(1),nn2(2)));
    
    % Update contour positions
    ssx = gamma*P(:,1) + Fext(:,1) ;% 
    ssy = gamma*P(:,2) + Fext(:,2) ;% 
    ssx(1) = ssx(1)+sf1(1); ssy(1) = ssy(1)+sf1(2);
    ssx(end) = ssx(end)+sf2(1); ssy(end) = ssy(end)+sf2(2);
    P(:,1) = B * ssx;
    P(:,2) = B * ssy;
    
    % Clamp contour to boundary
    P(:,1)=min(max(P(:,1),1),xdim);
    P(:,2)=min(max(P(:,2),1),ydim);
    
    Ps{ci}.pts = P;
end
  
P = cellInfoUpdate(Ps,I);