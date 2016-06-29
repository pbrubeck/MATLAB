function M = minus(M1,M2)
%MINUS  Subtract a scalar from a Moebius map.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: minus.m 298 2009-09-15 14:36:37Z driscoll $

M = M1 + (-M2);
