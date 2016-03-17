function [sem] = MathieuS(m, n, z, q)
j=mod(m,2);
k=(m+j)/2;
E(1:n-1)=q;
if j==0
    D=(2:2:2*n).^2;
else
    D=[1-q, (3:2:2*n-1).^2];  
end
[~,B]=trideigs(D,E);
B=B(:,k);
BB=zeros(2*n+2,1);
BB(3-j:2:end-j)=B;

plan=nfft(1,length(BB),numel(z));
plan.x=z(:)/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(BB);
nfft_trafo(plan);
sem=reshape(-imag(plan.f), size(z));
end