function R = polyMult(P, Q)
m=length(P);
n=length(Q);
R=zeros(1,m+n-1);
for i=1:m
    for j=1:n
        R(i+j-1)=R(i+j-1)+P(i)*Q(j);
    end
end
end
