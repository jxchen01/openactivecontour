function P=SnakeMoveIteration2D(B,P,Fext,gamma,kappa,delta,thickness)
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


% Update contour positions
ssx = gamma*P(:,1) + Fext1(:,1);% + Fext2(:,1);
ssy = gamma*P(:,2) + Fext1(:,2);% + Fext2(:,2);
P(:,1) = B * ssx;
P(:,2) = B * ssy;

% Clamp contour to boundary
P(:,1)=min(max(P(:,1),1),size(Fext,1));
P(:,2)=min(max(P(:,2),1),size(Fext,2));

    
