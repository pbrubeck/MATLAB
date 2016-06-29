function M = diskmap(M)
%DISKMAP Convert generic Schwarz-Christoffel map object to disk map.
%   DISKMAP(M) creates a diskmap object based on the polygon and
%   options contained in M.
%   
%   See the DISKMAP class documentation.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diskmap.m 298 2009-09-15 14:36:37Z driscoll $

M = diskmap(M.polygon,M.options);
