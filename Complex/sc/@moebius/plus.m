function M = plus(M1,M2)
%PLUS   Add a scalar to a Moebius map.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: plus.m 298 2009-09-15 14:36:37Z driscoll $

% Make the first one moebius.
if isa(M1,'double')
  tmp = M1;
  M1 = M2;
  M2 = tmp;
end

C = M1.coeff;
if isa(M2,'double') & length(M2)==1
  C(1) = C(1) + M2*C(3);
  C(2) = C(2) + M2*C(4);
else
  error('Addition not defined for these operands.')
end
M = moebius(C);
