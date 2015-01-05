function nP=headRegulation(P,I,flag,sz)

if(flag(1) && flag(2))
    tmp = P(2:end-1,:);
elseif(flag(1))
    tmp = P(2:end,:);
elseif(flag(2))
    tmp = P(1:end-1,:);
else
    error('error in head regulation');
end
nP=interpolateSingleContour(tmp,sz,size(P,1));