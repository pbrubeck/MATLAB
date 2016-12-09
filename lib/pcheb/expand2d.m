function UU = expand2d(U,Cx,Cy)
% EXPAND2D  Interpolate a matrix of collocation-point values.
%
% UU = EXPAND2D(U,C)
%    Interpolate a square matrix U of dimension n, using the specified
%    interpolant matrix C of size m-by-n, to yield UU of dimension m.
%
% UU = EXPAND2D(U,Cx,Cy)
%    Interpolate U (size ny-by-nx) using the specified interpolant
%    matrices Cx (mx-by-nx) and Cy (my-by-ny) to yield UU (my-by-mx).
%
% See also BARYCENTRIC, PBARY, PCHEB.

    if nargin < 3,  Cy = Cx;  end
    
    UU = Cy * U * Cx';

end
