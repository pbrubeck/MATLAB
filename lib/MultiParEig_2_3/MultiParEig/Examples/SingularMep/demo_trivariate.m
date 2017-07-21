%DEMO_TRIVARIATE   Demo for zeros of a system of trivariate polynomials 
%
% To solve a system of three cubic trivariate polynomials
%
% p(x,y,z) = 1 + 2x + 3y + 4z + 5x^2 + 6xy + 7xz + 8y^2 + 9yz + 10z^2 
%          + 11x^3 + 12x^2y + 13x^2z + 14xy^2 + 15xyz + 16xz^2 + 17y^3
%          + 18y^2z + 19yz^2 + 20z^3 = 0
% q(x,y,z) = 20 + 19x + 18y + 17z + 16x^2 + 15xy + 14xz + 13y^2 + 12yz 
%          + 11z^2 + 10x^3 + 9x^2y + 8x^2z + 7xy^2 + 6xyz + 5xz^2 + 4y^3 
%          + 3y^2z + 2yz^2 + z^3 = 0
% r(x,y,z) = 1 - 2x + 3y - 4z + 5x^2 - 6xy + 7xz - 8y^2 + 9yz - 10z^2 
%          + 11x^3 - 12x^2y + 13x^2z - 14xy^2 + 15xyz - 16xz^2 + 17y^3
%          - 18y^2z + 19yz^2 - 20z^3 = 0
%
% we linearize it as a singular three-parameter eigenvalue problem 
% with matrices of size 8 x 8 using a generalization of Linearization 1 
% from the reference. The finite eigenvalues are the zeroes. 
%
% See also: THREEPAREIG
%
% Reference: B. Plestenjak, M. E. Hochstenbach: Roots of bivariate polynomial 
% systems via determinantal representations, arXiv:1506.02291, 2015

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision 8.9.2015

% ------------------------------------------------------------------------
% p(x,y,z) = det(A1 + x*B1 + y*C1 + z*D1)
A1 = -eye(8);
A1(1,1) =  1;

B1 = zeros(8);
B1(1,:) = [2 5 0 0 11 14 16 0];
B1(2,1) = 1;  B1(5,2) = 1;  B1(8,3) = 1;

C1 = zeros(8);
C1(1,:) = [3 6 8 0 12 17 19 0];
C1(3,1) = 1;  C1(6,3) = 1; 

D1 = zeros(8);
D1(1,:) = [4 7 9 10 13 18 20 15];
D1(4,1) = 1;  D1(7,4) = 1;

% ------------------------------------------------------------------------
% q(x,y,z) = det(A2 + x*B2 + y*C2 + z*D2)
A2 = -eye(8);
A2(1,1) =  20;

B2 = zeros(8);
B2(1,:) = [19 16 0 0 10 7 5 0];
B2(2,1) = 1;  B2(5,2) = 1;  B2(8,3) = 1;

C2 = zeros(8);
C2(1,:) = [18 15 13 0 9 4 2 0];
C2(3,1) = 1;  C2(6,3) = 1; 

D2 = zeros(8);
D2(1,:) = [17 14 12 11 8 3 1 6];
D2(4,1) = 1;  D2(7,4) = 1;

% ------------------------------------------------------------------------
% r(x,y,z) = det(A3 + x*B3 + y*C3 + z*D3)
A3 = -eye(8);
A3(1,1) =  1;

B3 = zeros(8);
B3(1,:) = [-2 5 0 0 11 -14 -16 0];
B3(2,1) = 1;  B3(5,2) = 1;  B3(8,3) = 1;

C3 = zeros(8);
C3(1,:) = [3 -6 -8 0 -12 17 19 0];
C3(3,1) = 1;  C3(6,3) = 1; 

D3 = zeros(8);
D3(1,:) = [-4 7 9 -10 13 -18 -20 15];
D3(4,1) = 1;  D3(7,4) = 1;

% ------------------------------------------------------------------------
opts = [];
opts.singular = 1;
[x,y,z] = threepareig(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,opts);

zeroes = [x y z]