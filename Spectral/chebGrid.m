function x=chebGrid(n)
% Returns the extrema grid of the nth Chebyshev polynomial.
x=cos((0:n-1)'*pi/(n-1));
end