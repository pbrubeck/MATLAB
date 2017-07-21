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
% P. Holoborodko, Advanpix LLC.
% FreeBSD License, see LICENSE.txt

% PH 22.11.2016 : added precision-independency.
% PH 26.11.2016 : code simplifications and clean-ups.

% Last revision: 26.11.2016
narginchk(6, 7);

if nargin == 6 
    t = 1; 
end;

class_t = superiorfloat(A,B,C,D,E,F);

% Linearization
n = size(A,1);
Z = zeros(n,class_t);
I = eye(n,class_t);
P = [A       B   C;  Z I Z;  Z Z I];
Q = [Z       D t*E; -I Z Z;  Z Z Z];
R = [Z (1-t)*E   F;  Z Z Z; -I Z Z];
