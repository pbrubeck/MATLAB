function opt = scmapopt(M,varargin)
%SCMAPOPT Options structure for a Schwarz--Christoffel map object.
%   Same as the regular SCMAPOPT, but the first argument is the options
%   structure contained in the map M.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scmapopt.m 298 2009-09-15 14:36:37Z driscoll $

opt = scmapopt(M.options,varargin{:});
 