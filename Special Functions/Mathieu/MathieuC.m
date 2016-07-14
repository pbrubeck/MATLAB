function [cem] = MathieuC(m, q, z)
% Even Angular Mathieu Function
n=42;
A=MathieuA(m, q, n);
AA=zeros(2*n,1);
AA(mod(m,2)+1:2:end)=A;

plan=nfft(1,length(AA),numel(z));
plan.x=z(:)/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(AA);
nfft_trafo(plan);
cem=reshape(real(plan.f), size(z));
end