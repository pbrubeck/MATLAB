function [P1, P2] = blockprecond(cA,rA,m,n)
P1=zeros(m);
P2=zeros(n);
for i=1:m
    for j=1:n
        P1(:,i)=P1(:,i)+cA(i,j)/n; %Y(:,j)
        P2(:,j)=P2(:,j)+rA(i,j)/m; %Y(i,:)'
    end
end
end