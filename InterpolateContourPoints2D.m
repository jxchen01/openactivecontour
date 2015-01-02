function Ps=InterpolateContourPoints2D(Ps,nPoints,sz)
% This function resamples a few points describing a countour , to a smooth
% contour of uniformly sampled points.
%
% Ps=InterpolateContourPoints(Ps,nPoints)
%
% input,
%  Ps : Inpute Contours, N cells of size Ni x 2  
%  nPoints : Number of Contour points as output
% 
% output,
%  Ps : Uniformly sampled Contours, N cells of size nPoints x 2
%
% Function is written by D.Kroon University of Twente (July 2010)
% Modified by Jianxu Chen (University of Notre Dame) at Jan 2015


for i=1:1:numel(Ps)
    P=Ps{i}.pts;
    % Interpolate points inbetween
    r=min([4,floor((size(P,1)-1)/2)]);
    O(:,1)=interp(P(:,1),10,r);
    O(:,2)=interp(P(:,2),10,r);

    % Calculate distance between points
    dis=[0;cumsum(sqrt(sum((O(2:end,:)-O(1:end-1,:)).^2,2)))];

    % Resample to make uniform points
    K=zeros(nPoints,2);
    K(:,1) = interp1(dis,O(:,1),linspace(0,dis(end),nPoints));
    K(:,2) = interp1(dis,O(:,2),linspace(0,dis(end),nPoints));
    
    % Clamp contour to boundary
    K(:,1)=min(max(K(:,1),1),sz(1));
    K(:,2)=min(max(K(:,2),1),sz(2));

    Ps{i}.pts=K;
    
    clear O K dis r P
end


 