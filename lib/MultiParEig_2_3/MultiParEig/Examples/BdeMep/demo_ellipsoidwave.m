%DEMO_ELLIPSOIDWAVE  First eigenmodes of a tri-axial ellipsoid
%
% This example computes first 15 eigenmodes of an ellipsoid with semi-axes 1, 1.5, and 2
% The related Helmholtz equation separated in ellipsoidal coordinates is
% discretized with the Chebyshev collocation and solved as a
% three-parameter eigenvalue problem

% Reference: M. Willatzen and L. C. Lew Yan Voon, Numerical implementation of the ellipsoidal 
% wave equation and application to ellipsoidal quantum dots, Comput. Phys. Commun. 171 (2005) 1-18.
%
% See also: ELLIPSOIDWAVE_MEP, ELLIPSOID_EIGS, DEMO_ELLIPSOIDWAVE_FIGS.

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 21.11.2016 : correct rho sigma and tau output
% Last revision: 21.11.2016

N = 15; % number of collocation nodes 
x0 = 1;
y0 = 1.5;
z0 = 2;

a2 = z0^2-x0^2;
b2 = z0^2-y0^2;
c = a2/b2;
Momega = []; 
Mlambda = []; 
Mmu= []; 
Meta = [];
Mrst = [];

neig = 5;
RST = [0 0 0;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1; 1 1 1];
for izb = 1:8
    fprintf('Computing first %d eigenmodes of type (r,s,t) = (%d,%d,%d) \n',neig,RST(izb,1),RST(izb,2),RST(izb,3))
    [omega,X1,X2,X3,xi1,xi2,xi3,lambda,mu,eta] = ellipsoid_eigs(x0,y0,z0,RST(izb,1),RST(izb,2),RST(izb,3),neig,N,N,N);
    Momega = [Momega; omega];
    Mlambda = [Mlambda; lambda]; 
    Mmu = [Mmu; mu]; 
    Meta = [Meta; eta];
    Mrst = [Mrst; ones(neig,1)*RST(izb,:)];
end

[Momega, ord] = sort(Momega);
Mlambda = Mlambda(ord); 
Mmu = Mmu(ord);
Meta = Meta(ord); 
Mrst = Mrst(ord,:); 

disp('  omega      |   lambda       mu           eta        | rho sigma tau  ');
disp('-----------------------------------------------------------------------');
for j=1:15
   fprintf('%12.8f | %12.8f %12.8f %12.8f |  %d    %d    %d \n',...
       Momega(j), Mlambda(j),Mmu(j),Meta(j),Mrst(j,1),Mrst(j,2),Mrst(j,3))
end


