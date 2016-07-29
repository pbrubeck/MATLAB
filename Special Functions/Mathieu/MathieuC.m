function [cem] = MathieuC(m, q, z)
% Even Angular Mathieu Function ce_m(q; z)
if(~isreal(z))
    cem=MathieuJe(m, q, 1i*z);
    return;
end

n=42; s=mod(m,2);
A=MathieuA(m, q, n);
AA=zeros(4*n,1);
AA(1+s:2:2*n+s)=A;

plan=nfft(1,length(AA),numel(z));
plan.x=z(:)/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(AA);
nfft_trafo(plan);
cem=reshape(real(plan.f), size(z));
end