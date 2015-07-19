function m = Illinois(f, a, b)
% Solves for f(m)=0 over an initial interval [a,b].
% Modified false position.
side=0;
ym=1;
ya=f(a);
yb=f(b);
while(abs(ym)>eps) 
    m=(yb*a-ya*b)/(yb-ya);
    ym=f(m);
    if ym*ya<0
        b=m;
        yb=ym;
        if side==-1
            ya=ya/2;
        end
        side=-1;
    elseif ym*yb<0
        a=m;
        ya=ym;
        if side==1
            yb=yb/2;
        end
        side=1;
    else
        ym=0;
    end
end
end
