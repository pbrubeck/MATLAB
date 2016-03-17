function [cem] = MathieuC(m, n, z, q)
j=mod(m,2);
k=(m-j+2)/2;
E(1:n-1)=q;
if j==0
    D=(0:2:2*n-2).^2;
    E(1)=sqrt(2)*q;
else
    D=[1+q, (3:2:2*n-1).^2];  
end
[~,A]=trideigs(D,E);
A=A(:,k);
A=sign(A(k))*A;
A(1)=A(1)*(1-(1-j)/(2+sqrt(2)));
AA=zeros(2*n,1);
AA(1+j:2:end)=A;

plan=nfft(1,length(AA),numel(z));
plan.x=z(:)/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(AA);
nfft_trafo(plan);
cem=reshape(real(plan.f), size(z));
end