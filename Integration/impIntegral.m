function J = impIntegral(f, a, b)
% Allows evaluation of improper integrals using Romberg's method.
g=@(u) f(tan(u)).*(sec(u).^2);
J=Romberg(g, atan(a), atan(b));
end