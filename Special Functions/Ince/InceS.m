function [Sp] = InceS(p, m, q, z)
% Odd Ince Polynomial S^m_p(q; z)
B=InceB(p, m, q);
n=size(B,1); s=2-mod(p,2); k=2*(0:n-1)+s;
if isreal(z)
    B=padarray(B, max(0,5-length(B)), 'post');
    BB=[zeros(1+s,1); B(end:-1:2)/2; zeros(s,1); -B(1); -B(2:end)/2];
    plan=nfft(1,length(BB),numel(z));
    plan.x=z(:)/pi;
    nfft_precompute_psi(plan);
    plan.fhat=BB;
    nfft_trafo(plan);
    Sp=reshape(imag(exp(1i*s*z(:)).*plan.f), size(z));    
else
    Sp=reshape(sin(z(:)*k)*B, size(z));
end
end