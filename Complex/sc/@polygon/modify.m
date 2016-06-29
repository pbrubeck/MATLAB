function [p,indx] = modify(p)
%MODIFY Modify a polygon graphically.
%   See MODPOLY for usage instructions.

%   Copyright 1998 by Toby Driscoll.
%   $Id: modify.m 298 2009-09-15 14:36:37Z driscoll $

[w,beta,indx] = modpoly(vertex(p),angle(p)-1);
p = polygon(w,beta+1);
