function [P,Q,R] = linearize_quadtwopar(A,B,C,D,E,F,t)

%LINEARIZE_QUADTWOPAR  Linearize quadratic to linear two-parameter pencil
%
% [P,Q,R] = LINEARIZE_QUADTWOPAR(A,B,C,D,E,F,t)
% linearizes quadratic two-parameter pencil with matrices n x n
%   A + l*B + m*C + l^2*D + l*m*E + m^2*F
% as a linear two-parameter pencil with matrices 3n x 3n
%   P + l*Q + m*R 
% so that det(P + l*Q + m*R) = det(A + l*B + m*C + l^2*D + l*m*E + m^2*F)
% 
% real parameter t is a choice for the linearization (default is 1)
%
% See also: QUAD_TWOPAREIG.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision 8.9.2015

if nargin == 6 
    t = 1; 
end;

% Linearization
n = size(A,1);
Z = zeros(n);
I = eye(n);
P = [A       B   C;  Z I Z;  Z Z I];
Q = [Z       D t*E; -I Z Z;  Z Z Z];
R = [Z (1-t)*E   F;  Z Z Z; -I Z Z];
