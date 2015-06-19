function J = impIntegral(f, a, b)
a=atan(a);
b=atan(b);
g=@(u) f(tan(u))*(sec(u))^2;
J=Romberg(g, a, b);
end