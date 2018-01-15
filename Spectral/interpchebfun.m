function [ E ] = interpchebfun(n,xx)
% Interpolation operator from Chebyshev nodes to any other grid [-1,1]

% Non-equispaced FFT
xx=xx(:);
m=length(xx);
M=2*m-2;
N=2*n-2;
plan=nfft(1,N,M);
plan.x=acos([xx; -xx(end-1:-1:2)])/(2*pi);
nfft_precompute_psi(plan);

function [uhat]=adjoint(u)
    f=zeros(M,1); f(1:m)=u(:);
    plan.f=f;
    nfft_adjoint(plan);
    uhat=plan.fhat;
end

function [uu]=trafo(uhat)
    plan.fhat=uhat(:);
    nfft_trafo(plan);
    uu=plan.f(1:m);
end

function [uu]=Efun(u,tflag)
    if strcmp(tflag, 'transp')
        vhat=cellfun(@adjoint, num2cell(u, 1), 'UniformOutput', false);
        vv=fftshift(ifft((conj(cat(2, vhat{:}))),[],1));
        uu=vv(1:n,:); uu(2:end-1,:)=uu(2:end-1,:)+vv(N:-1:n+1,:);
        uu=real(uu);
    else
        dim=1;
        uhat=fftshift(ifft([u; u(end-1:-1:2,:)],[],dim,'symmetric'),dim); 
        uucell=cellfun(@trafo, num2cell(uhat, dim), 'UniformOutput', false);
        uu=ipermute(cat(ndims(u), uucell{:}), [dim, setdiff(1:ndims(u), dim)]);
        uu=real(uu);
    end
end
E=@Efun;
end