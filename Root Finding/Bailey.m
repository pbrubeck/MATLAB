function x = Bailey(f, ff, fff, x)
% Solves f(x)=0, requires f'(x), f''(x) and an initial aproximation.
y=f(x);
while(abs(y)>2*eps)
    yy=ff(x);
    x=x-y/(yy-y*fff(x)/(2*yy));
    y=f(x);
end
end