function Md = diff(M)
%DIFF   Differentiated SC map object.
%   DIFF(M) returns an object formally representing the derivative of
%   the map M.

%   Copyright 1998 by Toby Driscoll.
%   $Id: diff.m 298 2009-09-15 14:36:37Z driscoll $

Md = scmapdiff(M);
