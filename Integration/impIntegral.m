function J = impIntegral(f, a, b)
a=atan(a);
b=atan(b);
f=@(u) f(tan(u))*(sec(u))^2;
J=Romberg(f, a, b);
end