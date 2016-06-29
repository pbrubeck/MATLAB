function acc = accuracy(M)
%ACCURACY Apparent accuracy of Schwarz-Christoffel CR rectified map.

%   Copyright 1998 by Toby Driscoll.
%   $Id: accuracy.m 298 2009-09-15 14:36:37Z driscoll $

% Just return accuracy of underlying crdiskmap.
acc = accuracy(M.diskmap);
