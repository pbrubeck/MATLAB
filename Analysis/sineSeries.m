function b = sineSeries(f, P, n)
% Computes n coefficents of the cosine expansion of f over [a, b].
x=P*(0:n)/n;
h=f(x);
h(1)=(h(1)+h(end))/2;
b=2/n*dst(h(1:end-1));
end