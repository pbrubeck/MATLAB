function M = mtimes(M1,M2)
%MTIMES Multiply Moebius transformation by a scalar.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mtimes.m 298 2009-09-15 14:36:37Z driscoll $

% Make the first one moebius.
if isa(M1,'double')
  tmp = M1;
  M1 = M2;
  M2 = tmp;
end

C = M1.coeff;
if isa(M2,'double') & length(M2)==1
  C(1:2) = M2*C(1:2);
  M = moebius;
  M.coeff = C;
else
  error('Multiplication not defined for these operands.')
end
