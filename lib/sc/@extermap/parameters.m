function v = parameters(M)
%PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.

%   Copyright 1998 by Toby Driscoll.
%   $Id: parameters.m 298 2009-09-15 14:36:37Z driscoll $

v.prevertex = flipud(M.prevertex);
v.constant = M.constant;