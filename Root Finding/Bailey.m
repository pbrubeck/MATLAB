function x = Bailey(f, x)
% Solves f(x)=0, requires an initial aproximation.
x=ainit(x,2);
y=f(x);
while(abs(y{0})>10*eps)
    x=x-y{0}/(y{1}-y{0}*y{2}/(2*y{1}));
    y=f(x);
end
x=x{0};
end