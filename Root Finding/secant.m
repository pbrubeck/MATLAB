function b = secant(f, a, b)
% Solves f(x)=0 for x, requires two initial aproximations.
yb=1;
ya=f(a);
while(abs(yb)>eps)
    yb=f(b);
    temp=b;
    b=b-yb*(b-a)./(yb-ya);
    a=temp;
    ya=yb;
end
end
