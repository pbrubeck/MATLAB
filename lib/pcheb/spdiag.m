function A = spdiag(d)
%SPDIAG Sparse diagonal matrix.
%
% See also SPDIAGS, SPEYE, SPONE, SPARSE.

    d = d(:);
    n = length(d);
    A = spdiags(d,0,n,n);

end

