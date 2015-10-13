function [ out ] = RCcircuit(R, C, signal, t, mode)
% Returns the response of an RC to an input voltage signal
h=t(2)-t(1);
r=exp(-t/(R*C));
if(mode)
    s=signal;
    V=ifft(fft(s).*fft(r))*h; out=V;
else
    V=signal;
    s=ifft(fft(V)./fft(r))/h; out=s;
end
plot(t, [s; V]);
title(sprintf('R=%.2f, C=%.2f', R, C));
legend('Vi','Vo');
end