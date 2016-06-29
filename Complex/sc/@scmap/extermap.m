function M = extermap(M)
%EXTERMAP Convert generic Schwarz-Christoffel map object to exterior map.
%   EXTERMAP(M) creates a extermap object based on the polygon and
%   options contained in M.
%   
%   See the EXTERMAP class documentation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: extermap.m 298 2009-09-15 14:36:37Z driscoll $

M = extermap(M.polygon,M.options);
