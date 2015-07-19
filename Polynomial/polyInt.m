function J=polyInt(P, a, b)
n=size(P,2);
Q=zeros(size(P,1), n+1);
for i=1:n
    Q(:,i+1)=P(:,i)/i;
end
J=Horner(Q, b)-Horner(Q, a);
end