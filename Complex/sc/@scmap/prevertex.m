function z = prevertex(M)
%PREVERTEX Extract a vector of the prevertices of an S-C map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: prevertex.m 298 2009-09-15 14:36:37Z driscoll $

tmp = parameters(M);
if strmatch('prevertex',fieldnames(tmp))
  z = tmp.prevertex;
else
  msg = sprintf('Prevertices not defined for map of class %s\n',class(M));
  error(msg)
end
