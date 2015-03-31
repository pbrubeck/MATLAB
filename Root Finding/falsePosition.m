function m = falsePosition(f, a, b)
% Solves for f(m)=0 over an initial interval [a,b].

ii=0;
ym=0;
ya=f(a);
yb=f(b);
while(ii==0 || abs(ym)>1E-15) 
    m=a-ya*(b-a)/(yb-ya);
    ym=f(m);
    fprintf('i=%d \t a=%f \t b=%f \t m=%f \t f(m)=%f \n', ii, a, b, m, ym);
    if ym*ya<0
        b=m;
        yb=ym;
    elseif ym*yb<0
        a=m;
        ya=ym;
    else
        ym=0;
    end
    ii=ii+1;
end
