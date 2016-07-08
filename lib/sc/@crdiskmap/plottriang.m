function plottriang(M,varargin)
%PLOTTRIANG Plot triangulation for CR disk map.
%   PLOTTRIANG(M) plots the triangulation of the polygon used in CR disk
%   map M. PLOTTRIANG(M,LBL) also labels vertices and edges if LBL is
%   nonempty, 

%   Copyright 1998 by Toby Driscoll.
%   $Id: plottriang.m 298 2009-09-15 14:36:37Z driscoll $

plotptri(vertex(polygon(M)),M.qlgraph,varargin{:})
