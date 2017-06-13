function [ uu ] = interpcheb(xx,u)
N=size(u,1);
M=2*N-2;
uhat=ifft(u,M,'symmetric');

% Naive calculation
% uhat(1)=uhat(1)/2;
% uu=ChebT(2*uhat(1:N), xx);

% Non-equispaced FFT
plan=nfft(1,M,M);
plan.x=acos(xx([1:end,end-1:-1:2]))/(2*pi);
nfft_precompute_psi(plan);
plan.fhat=fftshift(uhat);
nfft_trafo(plan);
uu=plan.f;
uu=uu(1:N);
if isreal(u)
    uu=real(uu);
end
end