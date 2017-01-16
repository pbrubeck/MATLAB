function x = solveLaguerre(f, x)
% Solves f(x)=0, requires an initial aproximation.
x=ainit(x,2);
y=f(x);
while(abs(y{0})>10*eps)
    g=y{1}/y{0};
    h=g*g-y{2}/y{0};
    sgn=2*(g>=0)-1;
    a=2/(g+sgn*sqrt(2*h-g*g));
    x=x-a;
    y=f(x);
end
x=x{0};
end