function [Delta0,Delta1,Delta2] = twopar_delta(A1,B1,C1,A2,B2,C2)

%TWOPAR_DELTA   Delta matrices for a two-parameter eigenvalue problem
%
% [Delta0,Delta1,Delta2] = TWOPAR_DELTA(A1,B1,C1,A2,B2,C2) returns
% operator determinants 
%
% Delta0 = kron(C2, B1) - kron(B2, C1)
% Delta1 = kron(C2, A1) - kron(A2, C1)
% Delta2 = kron(A2, B1) - kron(B2, A1)
%
% related to the two-parameter eigenvalue problem
%
% A1 x = lambda B1 x + mu C1 x 
% A2 y = lambda B2 y + mu C2 y
%
% See also: THREEPAR_DELTA, MULTIPAR_DELTA.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

Delta0 = kron(B1,C2) - kron(C1,B2);
Delta1 = kron(A1,C2) - kron(C1,A2);
Delta2 = kron(B1,A2) - kron(A1,B2);