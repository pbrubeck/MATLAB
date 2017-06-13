function [uu] = interpcheb2(u,xx,yy)
% Interpolation from Chebyshev nodes to any other grid [-1,1]^2
N=size(u);
M=2*N(:)-2;
uhat=ifft(ifft(u,M(1),1,'symmetric'),M(2),2,'symmetric');

% Non-equispaced FFT
plan=nfft(2,M,prod(M));
xx=acos(xx([1:end, end-1:-1:2],[1:end, end-1:-1:2]))/(2*pi);
yy=acos(yy([1:end, end-1:-1:2],[1:end, end-1:-1:2]))/(2*pi);

plan.x=[xx(:), yy(:)];
nfft_precompute_psi(plan);
uhat=fftshift(uhat);
plan.fhat=uhat(:);
nfft_trafo(plan);

uu=reshape(plan.f, M(1), M(2));
uu=uu(1:N(1),1:N(2));
if isreal(u)
    uu=real(uu);
end
end

