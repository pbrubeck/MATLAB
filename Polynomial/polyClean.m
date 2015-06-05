function P=polyClean(P)
n=size(P,2);
if(n==1)
    return;
end
b=true;
while(n>1 && b)
    b=all(P(:,n)==0);
    n=n-1;
end
P=P(:,1:n+1);
end
