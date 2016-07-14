function [sem] = MathieuS(m, q, z)
% Odd Angular Mathieu Function
n=42;
B=MathieuB(m, q, 42);
r=mod(m,2);
BB=zeros(2*n+2,1);
BB(3-r:2:end-r)=B;

plan=nfft(1,length(BB),numel(z));
plan.x=z(:)/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(BB);
nfft_trafo(plan);
sem=reshape(-imag(plan.f), size(z));
end