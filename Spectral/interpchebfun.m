function [ E ] = interpchebfun(n,xx)
% Interpolation operator from Chebyshev nodes to any other grid [-1,1]

% Non-equispaced FFT
dim=1;
xx=xx(:);
m=length(xx);
N=2*n-2;
plan=nfct(1,n,m);
plan.x=acos(xx)/(2*pi);

function [uhat]=adjoint(u)
    plan.f=u(1:m);
    nfct_adjoint(plan);
    uhat=plan.fhat;
end

function [uu]=trafo(uhat)
    plan.fhat=uhat(1:n);
    nfct_trafo(plan);
    uu=plan.f;
end

function [uu]=Efun(u,tflag)
    if strcmp(tflag, 'transp')
        vhat=cellfun(@adjoint, num2cell(u,1), 'UniformOutput', false);
        vv=ifft(cat(2, vhat{:}),N,1,'symmetric');
        uu=vv(1:n,:); uu(2:end-1,:)=uu(2:end-1,:)+vv(end:-1:n+1,:);
        c=1+(m-1)*(xx(1)==1);
        uu(1,:)=uu(1,:)-u(c,:);
        uu(end,:)=uu(end,:)+u(c,:);
    else
        uhat=ifft(u,N,dim,'symmetric');
        uhat(2:n-1,:)=2*uhat(2:n-1,:);
        uucell=cellfun(@trafo, num2cell(uhat, dim), 'UniformOutput', false);
        uu=ipermute(cat(ndims(u), uucell{:}), [dim, setdiff(1:ndims(u), dim)]);
        uu(xx==-1,:)=u(end,:);
    end
end
E=@Efun;
end