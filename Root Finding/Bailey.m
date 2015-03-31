function xi = Bailey(f, ff, fff, xi)
% Solves f(x)=0, requires f'(x), f''(x) and an initial aproximation.
ii=0;
yi=0;
while(ii==0 || abs(yi)>1E-15)
    yi=f(xi);
    yyi=ff(xi);
    fprintf('i=%d \t xi=%f \t f(xi)=%f \n', ii, xi, yi);
    xi=xi-yi/(yyi-yi*fff(xi)/(2*yyi));
    ii=ii+1;
end
end
