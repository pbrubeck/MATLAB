%DEMO_ELLIPSOIDWAVE_JD   demo for the Jacobi-Davidson method for 
% three-parameter eigenvalue problems
%
% This example computes first 30 eigenmodes of an ellipsoid with semi-axes 1, 1.5, and 2
% The related Helmholtz equation separated in ellipsoidal coordinates is
% discretized with the Chebyshev collocation and solved as a
% three-parameter eigenvalue problem using the Jacobi-Davidson method

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

% Last revision: 21.11.2016

% We use matrices so large that we can not use threepareigs. We use threepareigs_jd 
% which calls threepareig to solve the smaller projected three-parameter.
% We compute twice as much eigenvalues with higher accuracy and faster as 
% in demo_ellipsoidwave.

N = 200; % number of collocation nodes 
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

neig = 10;
RST = [0 0 0;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1; 1 1 1];
for izb = 1:8
    fprintf('Computing first %d eigenmodes of type (r,s,t) = (%d,%d,%d) \n',neig,RST(izb,1),RST(izb,2),RST(izb,3))
    [omega,X1,X2,X3,xi1,xi2,xi3,lambda,mu,eta] = ellipsoid_eigs_jd(x0,y0,z0,RST(izb,1),RST(izb,2),RST(izb,3),neig,N,N,N);
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

disp('  omega          |   lambda           mu               eta            | rho sigma tau  ');
disp('-----------------------------------------------------------------------------------');
for j=1:30
   fprintf('%16.12f | %16.12f %16.12f %16.12f |  %d    %d    %d \n',...
       Momega(j), Mlambda(j),Mmu(j),Meta(j),Mrst(j,1),Mrst(j,2),Mrst(j,3))
end


