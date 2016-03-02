function [cem] = MathieuC(m, n, z, q)
j=mod(m,2);
E=q*ones(1,n-1);
if j==0
    D=(0:2:2*n-2).^2;
    E(1)=sqrt(2)*q;
else
    D=[1+q, (3:2:2*n-1).^2];  
end
[~,A]=trideigs(D,E);
A=A(:,(m+j)/2);
A=(2*(A(1)>0)-1)*A;
if j==0
    A=A/sqrt(A(1)^2+1);
end
k=numel(z);
AA=zeros(1,2*n-1);
AA(1+j:2:2*n-1+j)=A;
cem=k/2*ifft(AA,k,'symmetric')+(AA(1)-AA(end))/2;
end