function [P1, P2] = blockprecond(opA,m,n)
P1=zeros(m);
parfor i=1:m
    for j=1:n
        X=zeros(m,n);
        X(i,j)=1;
        Y=opA(X);
        P1(:,i)=P1(:,i)+Y(:,j)/n;
    end
end

P2=zeros(n);
parfor j=1:n
    for i=1:m
        X=zeros(m,n);
        X(i,j)=1;
        Y=opA(X);
        P2(:,j)=P2(:,j)+Y(i,:)'/m;
    end
end
end