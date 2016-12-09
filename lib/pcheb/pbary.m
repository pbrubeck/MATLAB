function [C,xx] = pbary(N,xel,xnh)
%PBARY  Barycentric interpolation for 1D p-element scheme.
%
% C = PBARY(N,xel,xx)
%   Interpolate from X=PCHEB(N,xel) to XX.
%   C is of size m-by-n, where n=length(x) and m=length(xx).
%  
% [C,xx] = PBARY(N,xel,xnh)
%   Uses vector XX returned by PCHEB(N,xel,xnh).
%
% See also PCHEB, BARYCENTRIC, PBARY_REMESH.

    [~,xx,C] = pcheb(N,xel,xnh);

end
