function [] = pulseCorr(n)
t=linspace(-2,3,n);
p=(t>=0)-(t>=1);
q=(t).*p;
r=sqrt(2*pi)/n*ifftshift(ifft(conj(fft(p)).*fft(q)));
plot(t,[p; q; r]);
end