function m = Illinois(f, a, b)
% Solves for f(m)=0 over an initial interval [a,b].
% Modified false position.
ii=0;
side=0;
ya=f(a);
yb=f(b);
ym=0;
while(ii==0 || abs(ym)>1E-15) 
    m=(yb*a-ya*b)/(yb-ya);
    ym=f(m);
    fprintf('i=%d \t a=%f \t b=%f \t m=%f \t y(m)=%f \n', ii, a, b, m, ym);
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
    ii=ii+1;
end
end
