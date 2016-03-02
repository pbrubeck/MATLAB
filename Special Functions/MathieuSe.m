function [sem] = MathieuSe(n, m, z, q)
j=mod(n,2);
E=q*ones(1,m-1);
if j==0
    D=(2:2:2*m).^2;
else
    D=[1-q, (3:2:2*m-1).^2];  
end
[~,B]=trideigs(D,E);
B=B(:,(n+j)/2);
%c=sqrt(m/2)*dst(B);
%c(1)=c(1)*sqrt(2);
sem=zeros(size(z));
for i=0:m-1
    sem=sem+B(i+1)*sin((2*i+2-j)*z);
end
end