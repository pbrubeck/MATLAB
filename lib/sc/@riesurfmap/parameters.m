function v = parameters(M)
%PARAMETERS Return a structure of the Schwarz-Christoffel map parameters.

%   Copyright 2002 by Toby Driscoll.
%   $Id: parameters.m 298 2009-09-15 14:36:37Z driscoll $

v.branch = M.branch;
v.prevertex = M.prevertex;
v.prebranch = M.prebranch;
v.constant = M.constant;
