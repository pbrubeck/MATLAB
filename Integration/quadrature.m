function J = quadrature(f, x, w)
% Computes the integral of a function given a quadrature rule.
J=f(x)*w(:);
end