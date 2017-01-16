function [A, G, H] = setBC(A, D, alpha, beta)
% Returns interior region differential operator suited for boundary 
% conditions of the form alpha(i)*u+beta(i)*u'=bc.
%
% The give-back matrix gives back the boundary values when multiplied by
% the interior values.
%
% u([1,N])=H*bc+G*u(2:N-1)
%
% References
% J. Hoepffner, Implementation of boundary conditions, 
% http://www.lmm.jussieu.fr/~hoepffner/boundarycondition.pdf

I=eye(size(A,1));
% Constraint matrix
C=diag(alpha)*I([1,end],:)+diag(beta)*D([1,end],:); 
% non-homogeneous contribution
H=inv(C(:,[1,end]));
% give-back matrix
G=-C(:,[1,end])\C(:,2:end-1); 
% Schur complement
A=A(2:end-1,2:end-1)+A(2:end-1,[1,end])*G;
end