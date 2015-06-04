function m = falsePosition(f, a, b)
% Solves for f(m)=0 over an initial interval [a,b].
i=0;
ym=1;
ya=f(a);
yb=f(b);
while(abs(ym)>eps) 
    m=a-ya*(b-a)/(yb-ya);
    ym=f(m);
    if ym*ya<0
        b=m;
        yb=ym;
    elseif ym*yb<0
        a=m;
        ya=ym;
    else
        ym=0;
    end
    i=i+1;
end
