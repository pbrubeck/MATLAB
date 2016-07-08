function [x,y] = vertex(p)
%VERTEX Vertices of a polygon.
%   VERTEX(P) returns the vertices of polygon P as a complex vector.
%   
%   [X,Y] = VERTEX(P) returns the vertices as two real vectors.
%   
%   See also POLYGON.

%   Copyright 1998 by Toby Driscoll.
%   $Id: vertex.m 298 2009-09-15 14:36:37Z driscoll $

x = p.vertex;
if nargout == 2
  y = imag(x);
  x = real(x);
end
