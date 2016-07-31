function [sem] = MathieuS(m, q, z)
% Odd Angular Mathieu Function se_m(q; z)
if q==0
    sem=sin(m*z);
    return;
end
n=42; s=2-mod(m,2); k=(s:2:2*n-2+s);
B=MathieuB(m, q, n);
if isreal(z)
    BB=zeros(4*n,1);
    BB(1+s:2:2*n+s)=B;
    plan=nfft(1,length(BB),numel(z));
    plan.x=z(:)/(2*pi);
    nfft_precompute_psi(plan);
    plan.fhat=fftshift(BB);
    nfft_trafo(plan);
    sem=reshape(-imag(plan.f), size(z));
else
    BB=B; BB(2:2:end)=-BB(2:2:end);
    b=sin(real(z)).^2>0.5;
    u=-2i*sqrt(q)*sin(z(b));
    w=2*sqrt(q)*cos(z(~b));
    sem=zeros(size(z));
    if s==2
        sem( b)=-(k*B)/B(1)/q*cot(z(b)).*sumj2k(B.*k', s, u);
        sem(~b)=(k*BB)/B(1)/q*tan(z(~b)).*sumj2k(BB.*k', s, w);
    else
        sem( b)=1i*(k*B)/B(1)/sqrt(q)*sumj2k(B, s, u);
        sem(~b)=sum(BB)/B(1)/sqrt(q)*tan(z(~b)).*sumj2k(BB.*k', s, w);
    end
end
end