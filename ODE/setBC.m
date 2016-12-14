function [A, G, CM] = setBC(A, D, alpha, beta)
% Returns interior region differential operator suited for homogenous
% boundary conditions of the form alpha(i)*u+beta(i)*u'=0.
%
% The give-back matrix gives back the boundary values when multiplied by
% the interior values.
%
% u([1,N])=G*u(2:N-1)
%
% References
% J. Hoepffner, Implementation of boundary conditions, 
% www.fukagata.mech.keio.ac.jp/~jerome/ (2007)

E=eye(size(A,1));

% Constraint Matrix
CM=diag(alpha)*E([1,end],:)+diag(beta)*D([1,end],:); 

% Give-back Matrix
G=-CM(:,[1,end])\CM(:,2:end-1); 

% Schur complement
A=A(2:end-1,2:end-1)+A(2:end-1,[1,end])*G;
end