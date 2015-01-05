function Ps=SnakeRegionUpdate(I,B,Ps,BMap,gamma,kappa,delta)
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
%     if(ci==5) % 4
%         keyboard;
%     end
    % retrieve info of current cell
    interiorIntensity = Ps{ci}.intensity;
    P =Ps{ci}.pts;
    R1=Ps{ci}.strip1;
    R2=Ps{ci}.strip2;
    NV=Ps{ci}.normvec;
    len=Ps{ci}.length;
    targetLength = Ps{ci}.targetLength;
    
    %%%%% prepare global image force %%%%%%%
    dphi_all = (I-interiorIntensity).^2-(I-exteriorIntensity).^2;
    %%%%% prepare global repelling force %%%%%%%
    id=setdiff(1:1:numContour+1,ci);
    rp_all=sum(RepelForce(:,:,id),3);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % forces exerted in normal direction of cell body
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % image force based on Chan-Vese model
    dphi(:,1) = interp2(dphi_all,R1(:,2),R1(:,1));
    dphi(:,2) = interp2(dphi_all,R2(:,2),R2(:,1));    
    dphi(isnan(dphi))=0;
    
    % repelling force
    rp(:,1) = interp2(rp_all,R1(:,2),R1(:,1));
    rp(:,2) = interp2(rp_all,R2(:,2),R1(:,1)); 
    rp(isnan(rp))=0;
    
    F_repel_Norm = -(rp(:,1)-rp(:,2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % forces exerted in tangential direction of head/tail
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nn=zeros(2,2);
    nn(1,:)=P(1,:)-P(2,:);%head
    nn(1,:)=nn(1,:)./hypot(nn(1,1),nn(1,2)); 
    nn(2,:)=P(end,:)-P(end-1,:); %tail
    nn(2,:)=nn(2,:)./hypot(nn(2,1),nn(2,2));
    
    %pp = P([1,end],:);
    pp=Ps{ci}.thickness.*nn + P([1,end],:);
    
    pp(:,1)=min(max(pp(:,1),1),xdim);
    pp(:,2)=min(max(pp(:,2),1),ydim);
    
    % repelling force 
    rp_head=interp2(rp_all,pp(:,2),pp(:,1));
    rp_head(isnan(rp_head))=0;
    
    % image force based on chan-vese model
    cv_head = interp2(dphi_all,pp(:,2),pp(:,1));
    cv_head(isnan(cv_head))=0;
    
    % the stretching force
    if(len<0.85*targetLength)
        phi=1;
    elseif(len<0.95*targetLength)
        phi=0.5;
    elseif(len>1.1*targetLength)
        phi=-1;
    else % expected length [0.95L, 1.05L]
        phi=0;
    end       
    sf=delta*phi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % normalize the image force along the contour
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cvMax=max([abs(dphi(:));abs(cv_head(:))]);
     
    %%%% normal direction %%%%
    
    % calculate the normalized image force in normal direction
    F_CV_Norm = -(dphi(:,1)-dphi(:,2))./cvMax;
    
    % Combing repeling force and image force into External Force Vector
    FN=F_CV_Norm+kappa.*F_repel_Norm;

    % upper bounded by +/- 1
    FN(FN>1)=1;
    FN(FN<-1)=-1;
    
    % make into vectors
    Fext(:,1)=NV(:,1).*FN ;
    Fext(:,2)=NV(:,2).*FN;

    %%%%%% tangential direction  %%%%%%%
    
    % calculate the normalized image force in tangential direction
    cv_head=cv_head./cvMax;
    
    % combine repelling force and image force
    ff=-(kappa.*rp_head+cv_head);
    
    % upper bounded by +/- 1
    ff(ff>1)=1; 
    ff(ff<-1)=-1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update contour positions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ss = gamma.*P + Fext ;
    ds = zeros(2,2);
    %ss(1,:) = ss(1,:) + ff(1).*nn(1,:);
    %ss(end,:) = ss(end,:) + ff(2).*nn(2,:);
    
    if(ff(1)>ff(2)) % larger means easiler to grow
        if(ff(1)<0)
            ds(1,:) = (ff(1)+2*sf).*nn(1,:);
        else
            ds(1,:) = (ff(1)+sf).*nn(1,:);
        end
        ds(2,:) = ff(2).*nn(2,:);
    else
        if(ff(2)<0)
            ds(2,:) = (ff(2)+2*sf).*nn(2,:);
        else
            ds(2,:) = (ff(2)+sf).*nn(2,:);
        end
        ds(1,:) = ff(1).*nn(1,:);
    end
    
    ss([1,end],:) = ss([1,end],:) + ds;

    P = B * ss;
    %P(:,2) = B * ssy;
    
    % Clamp contour to boundary
    P(:,1)=min(max(P(:,1),1),xdim);
    P(:,2)=min(max(P(:,2),1),ydim);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % regulation on head/tail
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flag=HighCurvature(P); % 1-by-2 logical array
    if(any(flag))
        P=headRegulation(P,I,flag,[xdim,ydim]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ps{ci}.pts = P;
end
