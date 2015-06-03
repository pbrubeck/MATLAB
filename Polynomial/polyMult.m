function R=polyMult(P, Q)
m=size(P,2);
n=size(Q,2);
R=zeros(max(size(P,1), size(Q,1)), m+n-1);
for i=1:m
    for j=1:n
        k=i+j-1;
        R(:,k)=R(:,k)+P(:,i).*Q(:,j);
    end
end
end