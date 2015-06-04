function xi = Bailey(f, ff, fff, xi)
% Solves f(x)=0, requires f'(x), f''(x) and an initial aproximation.
yi=1;
while(abs(yi)>eps)
    yi=f(xi);
    yyi=ff(xi);
    xi=xi-yi/(yyi-yi*fff(xi)/(2*yyi));
end
end
