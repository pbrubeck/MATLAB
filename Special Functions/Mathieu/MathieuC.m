function [cem] = MathieuC(m, q, z)
% Even Angular Mathieu Function ce_m(q; z)
if q==0
    cem=cos(m*z);
    return;
end
n=42; s=mod(m,2); k=(s:2:2*n-2+s);
A=MathieuA(m, q, n);
if isreal(z)
    AA=zeros(4*n,1);
    AA(1+s:2:2*n+s)=A;
    plan=nfft(1,length(AA),numel(z));
    plan.x=z(:)/(2*pi);
    nfft_precompute_psi(plan);
    plan.fhat=fftshift(AA);
    nfft_trafo(plan);
    cem=reshape(real(plan.f), size(z));
else
    AA=A; AA(2:2:end)=-AA(2:2:end);
    b=sin(real(z)).^2>0.5;
    u=-2i*sqrt(q)*sin(z(b));
    w=2*sqrt(q)*cos(z(~b));
    cem=zeros(size(z));
    if s==0
        cem( b)=sum(A)/A(1)*sumj2k(A, s, u);
        cem(~b)=sum(AA)/A(1)*sumj2k(AA, s, w);
    else
        cem( b)=1i*sum(A)/A(1)/sqrt(q)*cot(z(b)).*sumj2k(A.*k', s, u);
        cem(~b)=(k*AA)/A(1)/sqrt(q)*sumj2k(AA, s, w);
    end
end
end