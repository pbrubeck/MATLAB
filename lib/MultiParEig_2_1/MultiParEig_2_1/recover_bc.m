function U = recover_bc(X,G,k,r)

%RECOVER_BC  Recovers the removed degrees of freedom in BDE solution
%
% U = RECOVER_BC(X,G,k,r) reconstructs removed indices of a solution
% by applying give-back matrix and returns the complete solution (see
% Hoepffner)
%
% Input:
%   X : solution (on kept indices)
%   G : give-back matrix
%   k : kept indices
%   r : removed indiced
%
% See also: BDE2MEP, BDE3MEP

% Reference: J. Hoepffner, Implementation of boundary conditions, 
% www.fukagata.mech.keio.ac.jp/~jerome/ (2007)

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

N = length(k) + length(r); % number of variables
neig = size(X,2);  % number of eigenvectors

U = zeros(N,neig);
if ~isempty(r)
    U(r,:) = G*X;
end
U(k,:) = X;


