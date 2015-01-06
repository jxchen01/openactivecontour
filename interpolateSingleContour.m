function K=interpolateSingleContour(O, sz, nPoints)

% % Interpolate points inbetween
% r=min([4,floor((size(P,1)-1)/2)]);
% O(:,1)=interp(P(:,1),10,r);
% O(:,2)=interp(P(:,2),10,r);

% Calculate distance between points
dis=[0;cumsum(sqrt(sum((O(2:end,:)-O(1:end-1,:)).^2,2)))];

% Resample to make uniform points
K=zeros(nPoints,2);
K(:,1) = interp1(dis,O(:,1),linspace(0,dis(end),nPoints));
K(:,2) = interp1(dis,O(:,2),linspace(0,dis(end),nPoints));

% Clamp contour to boundary
K(:,1)=min(max(K(:,1),1),sz(1));
K(:,2)=min(max(K(:,2),1),sz(2));