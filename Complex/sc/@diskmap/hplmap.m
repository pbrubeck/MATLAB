function M1 = hplmap(M)
%HPLMAP Convert Schwarz-Christoffel disk map to a map from the half-plane.

%   Copyright 1998 by Toby Driscoll.
%   $Id: hplmap.m 298 2009-09-15 14:36:37Z driscoll $

p = polygon(M);
[z1,c1] = disk2hp(vertex(p),angle(p)-1,M.prevertex,M.constant);
M1 = hplmap(p,scmapopt(M),z1,c1);


