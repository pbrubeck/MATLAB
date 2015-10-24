function [pp] = noise(p,t)
p=p+randn(size(p));
p_hat=fft(p);
pp=sqrt(2*pi)/length(p)*ifft(conj(p_hat).*p_hat);
plot(t, [pp; cos(t)*sqrt(pi/2)]);
end

