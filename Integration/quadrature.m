function J = quadrature(f, x, w)
J=w*arrayfun(@(s) f(s), x)';
end
