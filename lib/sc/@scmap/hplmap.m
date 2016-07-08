function M = hplmap(M)
%HPLMAP Convert generic Schwarz-Christoffel map object to half-plane map.
%   HPLMAP(M) creates a hplmap object based on the polygon and
%   options contained in M.
%   
%   See the HPLMAP class documentation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hplmap.m 298 2009-09-15 14:36:37Z driscoll $

M = hplmap(M.polygon,M.options);
