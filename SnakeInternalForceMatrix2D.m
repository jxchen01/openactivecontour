function B=SnakeInternalForceMatrix2D(nPoints,alpha,beta,gamma)
%
% B=SnakeInternalForceMatrix2D(nPoints,alpha,beta,gamma)
%
% inputs,
%   nPoints : The number of snake contour points
%   alpha : membrame energy  (first order)
%   beta : thin plate energy (second order)
%   gamma : Step Size (Time)
%
% outputs,
%   B : The Snake Smoothness regulation matrix
%
% Function is written by D.Kroon University of Twente (July 2010)
beta=0;
% Penta diagonal matrix, one row:
b(1)=beta;
b(2)=-(alpha + 4*beta);
b(3)=(2*alpha + 6 *beta);
b(4)=b(2);
b(5)=b(1);

% Make the penta matrix (for every contour point)
A=b(1)*circshift(eye(nPoints),[2 0]);
A=A+b(2)*circshift(eye(nPoints),[1 0]);
A=A+b(3)*circshift(eye(nPoints),[0 0]);
A=A+b(4)*circshift(eye(nPoints),[-1 0]);
A=A+b(5)*circshift(eye(nPoints),[-2 0]);

% modification for open contour
% second row
A(2,1)=-alpha; A(2,2)=2*alpha; A(2,3)=-alpha; A(2,4)=0; A(2,end)=0;
% first row
A(1,1:3)=-0.5*A(2,1:3); A(1,end-1:end)=0;
% second last row
A(end-1,1)=0;A(end-1,end)=-alpha; A(end-1,end-1)=2*alpha; A(end-1,end-2)=-alpha; A(end-1,end-3)=0;
% last row
A(end,1:2)=0;A(end,end-2:end)=-0.5*A(end-1,end-2:end);

% Calculate the inverse
B=inv(A + gamma.* eye(nPoints));


