%DEMO_SINGULAR_MULTIPAREIG    Demo singular three-parameter eigenvalue problem with 2 x 2 matrices

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

A = cell(3,4);

% det(W1(x,y,z) = 1 + x + y + z;
A{1,1} = [1 0; 0 1]; 
A{1,2} = [1 0; 0 0]; 
A{1,3} = [1 0; 0 0]; 
A{1,4} = [1 0; 0 0];

% det(W2(x,y,z) = 1 + 2x + 3y + 4z;
A{2,1} = [1 0; 0 1]; 
A{2,2} = [2 0; 0 0]; 
A{2,3} = [3 0; 0 0]; 
A{2,4} = [4 0; 0 0];

% det(W3(x,y,z) = 1 + 2x - y + 2z;
A{3,1} = [1 0; 0 1]; 
A{3,2} = [2 0; 0 0]; 
A{3,3} = [-1 0; 0 0]; 
A{3,4} = [2 0; 0 0];

% The only regular solution is x = 4/3, y = 1/3, z = -2/3

opts=[];
opts.singular = 1;

[lambda,X,Y] = multipareig(A,opts);
lambda
