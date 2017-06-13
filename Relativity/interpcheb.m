function [ uu ] = interpcheb(u,xx,dim)
% Interpolation from Chebyshev nodes to any other grid [-1,1]
if nargin<3
    [~,dim]=max(size(u)>1);
end
xx=xx(:);
n=size(u,dim);
N=2*n-2;
uhat=fftshift(ifft(u,N,dim,'symmetric'),dim);

% Naive calculation
% uhat(1)=uhat(1)/2;
% uu=ChebT(2*uhat(1:n), xx);

% Non-equispaced FFT
m=length(xx);
M=2*m-2;
plan=nfft(1,N,M);
plan.x=acos(xx([1:end end-1:-1:2]))/(2*pi);
nfft_precompute_psi(plan);

order=[dim, 1:dim-1, dim+1:ndims(u)];
uhat=permute(uhat, order);
k=size(uhat);
uhat=reshape(uhat, N, []);
uu=zeros(m, size(uhat,2));
for i=1:size(uhat,2)
    plan.fhat=uhat(:,i);
    nfft_trafo(plan);
    if isreal(u)
        uu(:,i)=real(plan.f(1:m));
    else
        uu(:,i)=plan.f(1:m);
    end
end
uu=reshape(uu, [m, k(2:end)]);
uu=ipermute(uu, order);
end