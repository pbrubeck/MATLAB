function [cem] = MathieuC(m, n, z, q)
j=mod(m,2);
k=(m-j+2)/2;
E=q*ones(1,n-1);
if j==0
    D=(0:2:2*n-2).^2;
    E(1)=sqrt(2)*q;
else
    D=[1+q, (3:2:2*n-1).^2];  
end
[~,A]=trideigs(D,E);
A=A(:,k);
A=(2*(A(k)>0)-1)*A/sqrt(1+(1-j)*A(1)^2);
AA=zeros(1,2*n);
AA(1+j:2:end)=A;
cem=real(fft(AA,numel(z)));
end