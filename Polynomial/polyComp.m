function R=polyComp(P, Q)
m=size(P,2);
R=P(:,m);
for j=m-1:-1:1
    R=polyMult(R, Q);
    R(:,1)=R(:,1)+P(:,j);
end
R=polyClean(R);
end