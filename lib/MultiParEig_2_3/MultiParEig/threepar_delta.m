function [Delta0,Delta1,Delta2,Delta3] = threepar_delta(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3)

%THREEPAR_DELTA   Delta matrices for a three-parameter eigenvalue problem
%
% [Delta0,Delta1,Delta2,Delta3] = THREEPAR_DELTA(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3)
% returns operator determinants
% Delta0 = | B1 C1 D1; B2 C2 D2; B3 C3 D3 |
% Delta1 = | A1 C1 D1; A2 C2 D2; A3 C3 D3 |
% Delta2 = | B1 A1 D1; B2 A2 D2; B3 A3 D3 |
% Delta3 = | B1 C1 A1; B2 C2 A2; B3 C3 A3 |
% 
% related to the three-parameter eigenvalue problem
%
% A1 x1 = lambda B1 x1 + mu C1 x1 + eta D1 x1 
% A2 x2 = lambda B2 x2 + mu C2 x2 + eta D2 x2 
% A3 x3 = lambda B3 x3 + mu C3 x3 + eta D3 x3 
%
% See also: TWOPAR_DELTA, MULTIPAR_DELTA.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015
% BP 14.09.2016 : fixed bug in description

Delta0 = kron(B1,kron(C2,D3)-kron(D2,C3)) - kron(C1,kron(B2,D3)-kron(D2,B3)) + kron(D1,kron(B2,C3) - kron(C2,B3));
Delta1 = kron(A1,kron(C2,D3)-kron(D2,C3)) - kron(C1,kron(A2,D3)-kron(D2,A3)) + kron(D1,kron(A2,C3) - kron(C2,A3));
Delta2 = kron(B1,kron(A2,D3)-kron(D2,A3)) - kron(A1,kron(B2,D3)-kron(D2,B3)) + kron(D1,kron(B2,A3) - kron(A2,B3));
Delta3 = kron(B1,kron(C2,A3)-kron(A2,C3)) - kron(C1,kron(B2,A3)-kron(A2,B3)) + kron(A1,kron(B2,C3) - kron(C2,B3));
