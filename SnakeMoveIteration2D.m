function P=SnakeMoveIteration2D(B,P,Fext,gamma,kappa,delta,thickness,targetLength)
% This function will calculate one iteration of ribbon Snake movement
%
% P=SnakeMoveIteration2D(S,P,Fext,gamma,kappa)
%
% inputs,
%   B : Internal force (smoothness) matrix
%   P : The contour points N x 2;
%   Fext : External vector field (from image)
%   gamma : Time step
%   kappa : External (image) field weight
%   delta : stretching force due to length prior
%
% outputs,
%   P : The new (moved) contour points N x 2;
%
% Function is written by D.Kroon University of Twente (July 2010)

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));
P(:,2)=min(max(P(:,2),1),size(Fext,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get image force on the ribbon snake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get two outer strips of the ribbon 
[R1,R2]=getRibbon(P,thickness,size(Fext,1),size(Fext,2));

Fext1_R1(:,1)=kappa*interp2(Fext(:,:,1),R1(:,2),R1(:,1));
Fext1_R1(:,2)=kappa*interp2(Fext(:,:,2),R1(:,2),R1(:,1));
% Interp2, can give nan's if contour close to border
Fext1_R1(isnan(Fext1_R1))=0;

Fext1_R2(:,1)=kappa*interp2(Fext(:,:,1),R2(:,2),R2(:,1));
Fext1_R2(:,2)=kappa*interp2(Fext(:,:,2),R2(:,2),R2(:,1));
% Interp2, can give nan's if contour close to border
Fext1_R2(isnan(Fext1_R2))=0;

Fext1=Fext1_R1+Fext1_R2;

% Calculate the baloonforce on the contour points
%N=GetContourNormals2D(P);
%Fext2=delta*N;

% Calculate the stretching force
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
ssx = gamma*P(:,1) + Fext1(:,1);% + Fext2(:,1);
ssy = gamma*P(:,2) + Fext1(:,2);% + Fext2(:,2);
ssx(1) = ssx(1)+sf1(1); ssy(1) = ssy(1)+sf1(2);
ssx(end) = ssx(end)+sf2(1); ssy(end) = ssy(end)+sf2(2);
P(:,1) = B * ssx;
P(:,2) = B * ssy;

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));
P(:,2)=min(max(P(:,2),1),size(Fext,2));

    
