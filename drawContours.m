function h=drawContours(Ps,c,h,iter)

for i=1:1:numel(h)
    delete(h{i})
end

h=cell(1,numel(Ps));
for i=1:1:numel(Ps)
    P=Ps{i}.pts;
    h{i}=plot(P(:,2),P(:,1),'r.');
    plot(P(:,2),P(:,1),'-','Color',[c 1-c 0]);  
    title(['iteration: ',num2str(iter)]);
    drawnow
end