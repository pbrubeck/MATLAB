function a = cosineSeries(f, L, n)
% Computes n coefficents of the cosine expansion of f over [0, L].
x=L*(0:n)/n;
h=f(x);
h(1)=(h(1)+h(end))/2;
a=sqrt(2/n)*dct(h(1:end-1));
a(1)=a(1)/sqrt(2);
end