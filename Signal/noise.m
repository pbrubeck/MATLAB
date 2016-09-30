function [pp] = noise(p,t)
pr=p+randn(size(p));
p_hat=fft(pr);
pp=2/length(p)*ifft(conj(p_hat).*p_hat);
plot(t, [pp; p]);
end

