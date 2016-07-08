function M = mtimes(M,c)
%   Scale the image of a map by a complex constant.

% Usually this will be invoked by a child object, which needs to adjust
% its own constant as well.

%   Copyright (c) 1998 by Toby Driscoll.
%   $Id: mtimes.m 298 2009-09-15 14:36:37Z driscoll $

% May need to swap arguments
if isa(M,'double') & isa(c,'scmap')
  tmp = M;
  M = c;
  c = tmp;
end

M.polygon = c*M.polygon;
