function flag=HighCurvature(P)

flag=false(1,2);

a=P(1,:)-P(2,:);
b=P(3,:)-P(2,:);

if(dot(a,b)>0)
    flag(1)=true;
end

a=P(end,:)-P(end-1,:);
b=P(end-2,:)-P(end-1,:);

if(dot(a,b)>0)
    flag(2)=true;
end
