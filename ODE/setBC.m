function [A, G, H] = setBC(A, D, a, b, c)
% Returns interior region differential operator suited for boundary 
% conditions of the form a*u+b*u_x+c*u_t=bc.
%
% The give-back matrix gives back the boundary values when multiplied by
% the interior values.
%
% u([1,N])=H*bc+G*u(2:N-1)
%
% References
% J. Hoepffner, Implementation of boundary conditions, 
% http://www.lmm.jussieu.fr/~hoepffner/boundarycondition.pdf

if(nargin<5)
    c=0;
end

I=eye(size(A,1));
% Constraint matrix
C=diag(a)*I([1,end],:)+diag(b)*D([1,end],:)+diag(c)*A([1,end],:); 
% non-homogeneous contribution
H=inv(C(:,[1,end]));
% give-back matrix
G=-C(:,[1,end])\C(:,2:end-1); 
% Schur complement
A=A(2:end-1,2:end-1)+A(2:end-1,[1,end])*G;
end