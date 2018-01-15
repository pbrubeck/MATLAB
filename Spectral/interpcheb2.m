function [uu] = interpcheb2(u,xx,yy)
% Interpolation from Chebyshev nodes to any other grid [-1,1]^2
n=size(u);
N=2*n-2;
uhat=fftshift(ifft(ifft(u,N(1),1,'symmetric'),N(2),2,'symmetric'));

% Non-equispaced FFT
m=size(xx);
M=2*m-2;
plan=nfft(2,N(:),prod(M));
xx=acos(xx([1:end, end-1:-1:2],[1:end, end-1:-1:2]))/(2*pi);
yy=acos(yy([1:end, end-1:-1:2],[1:end, end-1:-1:2]))/(2*pi);

plan.x=[xx(:), yy(:)];
nfft_precompute_psi(plan);
plan.fhat=uhat(:);
nfft_trafo(plan);

uu=reshape(plan.f, M);
uu=uu(1:m(1),1:m(2));
if isreal(u)
    uu=real(uu);
end
end

