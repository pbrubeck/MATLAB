function [Cp] = InceC(p, m, q, z)
% Even Ince Polynomial C^m_p(q; z)
A=InceA(p, m, q);
n=size(A,1); s=mod(p,2); k=2*(0:n-1)+s;
if isreal(z)
    A=padarray(A, max(0,5-length(A)), 'post');
    AA=[zeros(1+s,1); A(end:-1:2)/2; zeros(s,1); A(1); A(2:end)/2];
    plan=nfft(1,length(AA),numel(z));
    plan.x=z(:)/pi;
    nfft_precompute_psi(plan);
    plan.fhat=AA;
    nfft_trafo(plan);
    Cp=reshape(real(exp(1i*s*z(:)).*plan.f), size(z));
else
    Cp=reshape(cos(z(:)*k)*A, size(z));
end
end