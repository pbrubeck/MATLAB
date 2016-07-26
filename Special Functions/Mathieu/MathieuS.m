function [sem] = MathieuS(m, q, z)
% Odd Angular Mathieu Function se_m(q; z)
if ~isreal(z)
    sem=-1i*MathieuJo(m, q, 1i*z);
    return;
end

n=42; s=2-mod(m,2);
B=MathieuB(m, q, 42);
BB=zeros(4*n,1);
BB(1+s:2:2*n+s)=B;

plan=nfft(1,length(BB),numel(z));
plan.x=z(:)/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(BB);
nfft_trafo(plan);
sem=reshape(-imag(plan.f), size(z));
end