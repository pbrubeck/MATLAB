function [ pp ] = autoCorr( p )
p_hat=fft(p);
pp=sqrt(2*pi)/length(p)*ifftshift(ifft(conj(p_hat).*p_hat));
end