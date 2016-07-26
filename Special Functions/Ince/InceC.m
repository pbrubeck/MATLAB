function [Cp] = InceC(p, m, q, z)
% Even Ince Polynomial C^m_p(q; z)
A=InceA(p, m, q);
n=size(A,1); s=mod(p,2);
AA=zeros(max(8, 4*n), 1);
AA(1+s:2:2*n+s)=A;

plan=nfft(1,length(AA),numel(z));
plan.x=z(:)/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(AA);
nfft_trafo(plan);
Cp=reshape(real(plan.f), size(z));
end