function J = quadrature(f, x, w)
J=arrayfun(@(s) f(s), x)*w';
end