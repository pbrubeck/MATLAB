function [y0,y1,y2,y3] = threepar_delta_mv(z,A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3)

%THREEPAR_DELTA_MV   Efficient multiplication with Delta matrices for 3p
%
% [y0,y1,y2,y3] = threepar_delta_mv(z,A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3)
% returns products yk = Deltak * z for k=0,1,2,3, where operator 
% determinants
% Delta0 = | B1 C1 D1; B2 C2 D2; B3 C3 D3 |
% Delta1 = | A1 C1 D1; A2 C2 D2; A3 C3 D3 |
% Delta2 = | B1 A1 D1; B2 A2 D2; B3 A3 D3 |
% Delta3 = | B1 C1 A1; B2 C2 A2; B3 C3 A3 |
%
% are related to the three-parameter eigenvalue problem
%
% A1 x1 = lambda B1 x1 + mu C1 x1 + eta D1 x1 
% A2 x2 = lambda B2 x2 + mu C2 x2 + eta D2 x2 
% A3 x3 = lambda B3 x3 + mu C3 x3 + eta D3 x3 
 
% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 14.09.2016

n1 = size(A1,1);
n2 = size(A2,1);
n3 = size(A3,1);
 
y0 = mode123(z,B1,C2,D3) + mode123(z,C1,D2,B3) + mode123(z,D1,B2,C3) - mode123(z,D1,C2,B3) - mode123(z,C1,B2,D3) - mode123(z,B1,D2,C3); 
y1 = mode123(z,A1,C2,D3) + mode123(z,C1,D2,A3) + mode123(z,D1,A2,C3) - mode123(z,D1,C2,A3) - mode123(z,C1,A2,D3) - mode123(z,A1,D2,C3); 
y2 = mode123(z,B1,A2,D3) + mode123(z,A1,D2,B3) + mode123(z,D1,B2,A3) - mode123(z,D1,A2,B3) - mode123(z,A1,B2,D3) - mode123(z,B1,D2,A3); 
y3 = mode123(z,B1,C2,A3) + mode123(z,C1,A2,B3) + mode123(z,A1,B2,C3) - mode123(z,A1,C2,B3) - mode123(z,C1,B2,A3) - mode123(z,B1,A2,C3); 
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary function for mode-1 matrix - tensor multipliation  
function w = mode1(z,A,n1,n2,n3) 
    w = reshape((A*reshape(z,[n2*n3 n1]).').',[n1*n2*n3 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary function for mode-2 matrix - tensor multipliation  
function w = mode2(z,A,n1,n2,n3) 
    w = reshape(reshape(A*reshape(reshape(z,[n3 n1*n2]).',[n2 n1*n3]),[n1*n2 n3]).',[n1*n2*n3 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary function for mode-3 matrix - tensor multipliation  
function w = mode3(z,A,n1,n2,n3) 
    w = reshape(A*reshape(z,[n3 n1*n2]),[n1*n2*n3 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary function for multiplication (A(x)B(x)C)*z
function w = mode123(z,A,B,C)
    n1 = size(A,1);
    n2 = size(B,1);
    n3 = size(C,1);
    w = mode1(z,A,n1,n2,n3);
    w = mode2(w,B,n1,n2,n3);
    w = mode3(w,C,n1,n2,n3);
end