function M = uminus(M)
%UMINUS Negate a Moebius transformation.
%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: uminus.m 298 2009-09-15 14:36:37Z driscoll $

M.coeff(1:2) = -M.coeff(1:2);
