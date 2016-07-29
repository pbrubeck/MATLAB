function [Sp] = InceS(p, m, q, z)
% Odd Ince Polynomial S^m_p(q; z)
B=InceB(p, m, q);
n=size(B,1); s=2-mod(p,2);
if(~isreal(z))
    Sp=reshape(sin(z(:)*(2*(0:n-1)+s))*B, size(z));
    return;
end

BB=zeros(max(16, 4*n), 1);
BB(1+s:2:2*n+s)=B;
plan=nfft(1,length(BB),numel(z));
plan.x=z(:)/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(BB);
nfft_trafo(plan);
Sp=reshape(-imag(plan.f), size(z));
end