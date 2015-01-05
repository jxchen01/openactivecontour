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

    Ps{i}.pts= interpolateSingleContour(Ps{i}.pts,sz, nPoints);
    
end


 