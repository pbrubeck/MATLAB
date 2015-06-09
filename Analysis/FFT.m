function [X] = FFT(x, m, N, s)
%FFT Summary of this function goes here
%   Detailed explanation goes here
if(nargin==1)
    m=1;
    N=length(x);
    s=1;
end

if(N==1)
    r=length(x);
    X=x(1+mod(r-m, r));
    return;
end

E=FFT(x, m, N/2, 2*s);
O=FFT(x, m+s, N/2, 2*s);
X=zeros(1,N);
wn=exp(2i*pi/N);
tw=1;

for k=1:N/2
    t=tw*O(k);
    X(k)=E(k)+t;
    X(k+N/2)=E(k)-t;
    tw=tw*wn;
end

end