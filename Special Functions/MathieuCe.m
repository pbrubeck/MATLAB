function [cem] = MathieuCe(n, m, z, q)
j=mod(n,2);
E=q*ones(1,m-1);
if j==0
    D=(0:2:2*m-2).^2;
    E(1)=sqrt(2)*q;
else
    D=[1+q, (3:2:2*m-1).^2];  
end
[~,A]=trideigs(D,E);
A=A(:,(n+j)/2);
if j==0
    A=A/sqrt(A(1)^2+1);
end
%c=sqrt(m/2)*dct(A);
%c(1)=c(1)*sqrt(2);
cem=zeros(size(z));
for i=0:m-1
    cem=cem+A(i+1)*cos((2*i+j)*z);
end
end