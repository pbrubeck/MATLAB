function fp = evaldiff(M,zp)
%EVALDIFF Derivative of Schwarz-Christoffel disk map at points.
%   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
%   disk map M at the points ZP.
%   
%   See also DISKMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: evaldiff.m 298 2009-09-15 14:36:37Z driscoll $

z = M.prevertex;
c = M.constant;
beta = angle(polygon(M)) - 1;

fp = dderiv(zp,z,beta,c);
