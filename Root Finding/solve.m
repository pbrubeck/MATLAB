function x = solve(f, ff, fff, x)
% Solves f(x)=0, requires f'(x), f''(x) and an initial aproximation.
y=f(x);
while(abs(y)>2*eps)
    G=ff(x)/y;
    H=G*G-fff(x)/y;
    sgn=2*(G>=0)-1;
    a=2/(G+sgn*sqrt(2*H-G*G));
    x=x-a;
    y=f(x);
end
end