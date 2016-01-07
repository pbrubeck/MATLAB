function [b] = fpm(a, Q)
% Fast Polynomial Multiplication R(x) = P(x)*Q(x)
% P(x) of degree n, given as a Chebyshev series a[k]*T[k](x)
% Q(x) of degree m, must be sampled on cos(pi*(1:2:2*M-1)/(2*M))
% where M/2 <= m+n < M = 2^s, i.e. M := 2^floor(1+log2(m+n))
% R(x) is returned as a Chebyshev series b[k]*T[k](x)
a(1)=a(1)*sqrt(2);
b=dct(idct(a, length(Q)).*Q); 
b(1)=b(1)/sqrt(2);
end