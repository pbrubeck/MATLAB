function fp = evaldiff(M,zp)
%EVALDIFF Derivative of Schwarz-Christoffel rectangle map at points.
%   EVALDIFF(M,ZP) computes the derivative of the Schwarz-Christoffel
%   rectangle map M at the points ZP.
%   
%   See also RECTMAP, EVAL.

%   Copyright 1998 by Toby Driscoll.
%   $Id: evaldiff.m 298 2009-09-15 14:36:37Z driscoll $

z = M.prevertex;
c = M.constant;
beta = angle(polygon(M)) - 1;

fp = rderiv(zp,z,beta,M.constant,M.stripL);
