function M = stripmap(M)
%STRIPMAP Convert generic Schwarz-Christoffel map object to strip map.
%   STRIPMAP(M) creates a stripmap object based on the polygon and
%   options contained in M.
%   
%   See the STRIPMAP class documentation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stripmap.m 298 2009-09-15 14:36:37Z driscoll $

M = stripmap(M.polygon,M.options);
