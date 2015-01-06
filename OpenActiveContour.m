function [P,J]=OpenActiveContour(I,P,Options)

% [O,J]=OpenActiveContour(I,P,Options)
%  
% inputs,
%   I : An Image of type double preferable ranged [0..1]
%   P : List with coordinates descriping the rough contour N x 2
%   Options : A struct with all snake options
%   
% outputs,
%   O : List with coordinates of the final contour M x 2
%   J : Binary image with the segmented region
%
% options (general),
%  Option.Verbose : If true show important images, default false
%  Options.nPoints : Number of contour points, default 100
%  Options.Gamma : Time step, default 1
%  Options.Iterations : Number of iterations, default 100
%
% options (internal force)
%  Options.Alpha : Membrame energy  (first order), default 0.2
%  Options.Beta : Thin plate energy (second order), default 0.0

% options (Snake)
%  Options.Delta : stretching force due to length prior
%  Options.Kappa : Weight of repelling force, default 0.2

% Function is written by D.Kroon University of Twente (July 2010)
% Modified by Jianxu Chen (University of Notre Dame) at Jan 2015

% Process inputs
defaultoptions=struct('Verbose',false,'nPoints',25,'Alpha',0.2,'Beta',0.0,'Delta',2,...
    'Gamma',1,'Kappa',0.2,'Iterations',100);

if(~exist('Options','var')), 
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))), 
        warning('snake:unknownoption','unknown options found');
    end
end

% If color image convert to grayscale
if(size(I,3)==3)
    I=rgb2gray(I); 
else
    I=mat2gray(I);
end

% Make an uniform sampled contour description
P=InterpolateContourPoints2D(P,Options.nPoints,size(I));
P=cellInfoUpdate(P,I);
if(Options.Verbose)
    figure(2), imshow(I), hold on; myHandle=drawContours(P,0,[],0);
end

% Make the interal force matrix (smooth the contour)
S=SnakeInternalForceMatrix2D(Options.nPoints,Options.Alpha,Options.Beta,Options.Gamma);

% Transform the Image into an External Energy Image
Eext = ExternalForceImage2D(I,0.04, 2, 0.01 ,6);
 
% Make the external force (flow) field.
Fx=ImageDerivatives2D(Eext,20,'x');
Fy=ImageDerivatives2D(Eext,20,'y');
Fext(:,:,1)=-Fx*2*20^2;
Fext(:,:,2)=-Fy*2*20^2;

Fext=GVFOptimizeImageForces2D(Fext, 0.1, 5, 1.0);

for i=1:Options.Iterations    
    %P=SnakeMoveIteration2D(S,P,Fext,Options.Gamma,Options.Kappa,Options.Delta,thickness,targetLength);
    P=SnakeRegionUpdate(I,S,P,Fext,zeros(size(I)),Options.Gamma,Options.Kappa,Options.Delta);
    
    %if(mod(i,5)==0)
        P=InterpolateContourPoints2D(P,Options.nPoints,size(I));
    %end
    P = cellInfoUpdate(P,I);
    
    % Show current contour
    if(Options.Verbose)
        myHandle=drawContours(P,i/Options.Iterations,myHandle,i);
    else
        disp(['iteration: ',num2str(i)]);
    end
    
    if(stopCheck(P))
        break;
    end
end

if(nargout>1)
     J=DrawSegmentedArea2D(P,I);
end

end

function flag=stopCheck(P)
    flag=true;
    for i=1:1:numel(P)
        %disp([abs(P{i}.length - P{i}.targetLength),0.05*P{i}.targetLength]);
        if(abs(P{i}.length - P{i}.targetLength)> 0.05*P{i}.targetLength)
            flag=false;
            break;
        end
    end
end

