function [ uu ] = interpnfft(x,u,xx,dim)
% Interpolation from u(x) to uu(xx)
if nargin<4
    [~,dim]=max(size(u)>1);
end

% Non-equispaced FFT
x=x(:);   n=length(x);  N=2*n-2;
xx=xx(:); m=length(xx); M=2*m-2;

% Retrieve coefficients: adjoint transform N -> N
plan1=nfft(1,N,N);
plan1.x=[acos(x(1:end)); -acos(x(end-1:-1:2))]/(2*pi);
nfft_precompute_psi(plan1);

% Evaluation: direct transform N -> M
plan2=nfft(1,N,M);
plan2.x=[acos(xx(1:end)); -acos(xx(end-1:-1:2))]/(2*pi);
nfft_precompute_psi(plan2);

function [uu]=dointerp(u)
    u=u(:);
    plan1.f=u([1:end end-1:-1:2]);
    nfft_adjoint(plan1);
    plan2.fhat=plan1.fhat/N;
    nfft_trafo(plan2);
    uu=plan2.f(1:m);
end

uucell=cellfun(@dointerp, num2cell(u, dim), 'UniformOutput', false);
uu=ipermute(cat(ndims(u), uucell{:}), [dim, setdiff(1:ndims(u), dim)]);

if isreal(u)
    uu=real(uu);
end
end