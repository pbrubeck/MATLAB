function h = plotcdt(p,T,varargin)
%PLOTCDT Plot constrained Delaunay triangulation.
%   PLOTCDT(P,T) plots the CDT of P computed by CDT. PLOTCDT(P,T,1) labels
%   the edges and vertices.
%   
%   H = PLOTCDT(P,T) returns a vector of handles for the edges.
%   
%   See also CDT.

%   Copyright 1998 by Toby Driscoll.
%   $Id: plotcdt.m 298 2009-09-15 14:36:37Z driscoll $

han = plotptri(p.vertex,T.edge,varargin{:});

if nargout > 0
  h = han;
end
