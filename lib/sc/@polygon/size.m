function n = size(p,m)
%   Number of vertices.

%   Copyright 1998 by Toby Driscoll.
%   $Id: size.m 298 2009-09-15 14:36:37Z driscoll $

if nargin == 1
  n = [length(p.vertex) 1];
elseif m == 1
  n = length(p.vertex);
else
  n = 1;
end

  