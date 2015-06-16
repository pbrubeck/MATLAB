function J = quadrature(f, x, w)
% Computes the integral of a function given a quadrature rule.
J=arrayfun(@(s) f(s), x)*w(:);
end