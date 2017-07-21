%DEMO_ELLIPSOIDWAVE_CONVERGENCE   convergence of approximations for 
% ellipsoid eigenmodes by Chebyshev collocation 
%
% This example computes first 20 eigenmodes of an ellipsoid with semi-axes 1, 1.5, and 2
% The related Helmholtz equation separated in ellipsoidal coordinates is
% discretized with the Chebyshev collocation and solved as a
% three-parameter eigenvalue problem using the Jacobi-Davidson method. The
% result is refined in multiple precision. Eigenmodes for different number
% of Chebyshev points are computed: 10, 15, 20, 25, 30, 100, differences
% that show the convergence are presented in the table. We also show the
% errors of approximations obtained in double precision.
%
% This example requires Multiprecision Computing Toolbox for MATLAB, see
% http://www.advanpix.com/
%
% See also: ELLIPSOIDWAVE_MEP, ELLIPSOID_EIGS, DEMO_ELLIPSOIDWAVE_FIGS.

% References: 
%  - M. Willatzen and L. C. Lew Yan Voon, Numerical implementation of the ellipsoidal 
%    wave equation and application to ellipsoidal quantum dots, Comput. Phys. Commun. 171 (2005) 1-18.
%  - B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
%    collocation for multiparameter eigenvalue problems arising from separable
%    boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 01.12.2016

% Example output:
% 
%   "Exact" eigenmodes and errors for "exact" approximations obtained by Chebyshev collocation
% ------------------------------------------------------------------------------------------------------
%   omega (N=100)                | error(N=30) | error(N=25) | error(N=20) | error(N=15) | rho sigma tau  
% ------------------------------------------------------------------------------------------------------
%   2.34458979280883788429599017 |    1.23e-80 |    9.91e-64 |    1.19e-47 |    1.34e-32 |  0    0    0 
%   2.94367434748453034138116762 |    1.02e-75 |    9.93e-60 |    1.48e-44 |    2.19e-30 |  1    0    0 
%   3.20795093379006047099453100 |    1.63e-73 |    6.73e-58 |    4.27e-43 |    2.69e-29 |  0    1    0 
%   3.57728276505342269224500042 |    1.24e-69 |    1.40e-54 |    2.30e-40 |    3.41e-27 |  0    0    0 
%   3.78641651063680513812657396 |    4.74e-70 |    4.36e-55 |    6.31e-41 |    9.47e-28 |  1    1    0 
%   3.82663625847381319117799359 |    7.96e-69 |    5.68e-54 |    6.25e-40 |    6.88e-27 |  0    0    1 
%   4.13064732283380355067391621 |    8.63e-66 |    2.33e-51 |    9.15e-38 |    3.28e-25 |  0    0    0 
%   4.23215870790611470111743051 |    3.90e-66 |    9.53e-52 |    3.48e-38 |    1.20e-25 |  1    0    0 
%   4.38693776442480328648309932 |    4.40e-66 |    9.32e-52 |    3.12e-38 |    1.09e-25 |  1    0    1 
%   4.38859178085098994822348040 |    2.94e-65 |    5.01e-51 |    1.27e-37 |    3.04e-25 |  0    1    0 
%   4.61577934035241530913177864 |    8.35e-65 |    1.07e-50 |    2.17e-37 |    4.61e-25 |  0    1    1 
%   4.70777812027019249248103697 |    3.13e-63 |    2.65e-49 |    3.36e-36 |    4.04e-24 |  1    0    0 
%   4.89789930584256271603929900 |    1.63e-61 |    7.38e-48 |    4.69e-35 |    2.50e-23 |  0    0    0 
%   4.97229440613881937292042883 |    5.81e-62 |    2.85e-48 |    2.08e-35 |    1.43e-23 |  0    0    1 
%   5.00681461492429543112942680 |    1.34e-62 |    7.04e-49 |    5.67e-36 |    4.46e-24 |  1    1    0 
%   5.07658292876947251862780778 |    2.65e-61 |    1.06e-47 |    6.36e-35 |    3.64e-23 |  0    1    0 
%   5.16933801956031837479064992 |    1.39e-62 |    6.63e-49 |    5.14e-36 |    4.33e-24 |  1    1    1 
%   5.30934082829887973029434192 |    2.40e-59 |    4.89e-46 |    1.39e-33 |    3.33e-22 |  0    0    0 
%   5.37807222799613368985795503 |    1.21e-58 |    2.35e-45 |    6.66e-33 |    1.74e-21 |  0    0    0 
%   5.44699842002785581489333857 |    1.77e-59 |    3.52e-46 |    1.05e-33 |    3.00e-22 |  0    0    1 
% ------------------------------------------------------------------------------------------------------
%  
%   "Exact" eigenmodes and errors for approximations obtained by Chebyshev collocation in double precision
% ---------------------------------------------------------------------------------------------------------------------
%   omega (N=100)                | error(N=100) | error(N=30) | error(N=25) | error(N=20) | error(N=15) | rho sigma tau  
% ---------------------------------------------------------------------------------------------------------------------
%   2.34458979280883788429599017 |     1.10e-12 |    8.12e-12 |    2.17e-12 |    4.78e-14 |    2.23e-12 |  0    0    0 
%   2.94367434748453034138116762 |     1.71e-09 |    1.13e-11 |    3.23e-12 |    9.50e-15 |    2.65e-16 |  1    0    0 
%   3.20795093379006047099453100 |     2.80e-10 |    1.09e-11 |    2.32e-11 |    7.70e-12 |    8.24e-14 |  0    1    0 
%   3.57728276505342269224500042 |     4.19e-12 |    7.55e-14 |    1.81e-13 |    6.48e-14 |    5.77e-15 |  0    0    0 
%   3.78641651063680513812657396 |     2.50e-12 |    4.46e-13 |    2.18e-12 |    4.06e-13 |    1.21e-14 |  1    1    0 
%   3.82663625847381319117799359 |     7.06e-10 |    3.91e-14 |    6.40e-12 |    2.50e-13 |    2.08e-14 |  0    0    1 
%   4.13064732283380355067391621 |     5.80e-11 |    6.01e-13 |    5.74e-14 |    6.25e-14 |    4.76e-13 |  0    0    0 
%   4.23215870790611470111743051 |     2.46e-12 |    1.47e-14 |    3.11e-15 |    4.84e-14 |    2.62e-14 |  1    0    0 
%   4.38693776442480328648309932 |     9.81e-10 |    1.66e-11 |    4.87e-14 |    4.28e-14 |    1.08e-15 |  1    0    1 
%   4.38859178085098994822348040 |     2.48e-09 |    1.26e-12 |    9.95e-14 |    4.70e-14 |    4.79e-14 |  0    1    0 
%   4.61577934035241530913177864 |     6.07e-11 |    3.05e-12 |    2.25e-12 |    2.86e-14 |    1.88e-14 |  0    1    1 
%   4.70777812027019249248103697 |     2.98e-11 |    1.05e-11 |    2.60e-14 |    2.02e-14 |    1.49e-14 |  1    0    0 
%   4.89789930584256271603929900 |     2.28e-09 |    2.38e-11 |    2.48e-14 |    6.85e-14 |    7.56e-14 |  0    0    0 
%   4.97229440613881937292042883 |     1.51e-10 |    4.78e-14 |    4.54e-14 |    5.61e-14 |    2.21e-14 |  0    0    1 
%   5.00681461492429543112942680 |     1.04e-09 |    9.64e-14 |    2.63e-14 |    3.60e-14 |    3.92e-15 |  1    1    0 
%   5.07658292876947251862780778 |     3.19e-13 |    3.86e-14 |    6.73e-15 |    1.99e-14 |    1.43e-13 |  0    1    0 
%   5.16933801956031837479064992 |     2.53e-10 |    1.98e-12 |    9.78e-16 |    1.87e-15 |    2.83e-14 |  1    1    1 
%   5.30934082829887973029434192 |     1.24e-09 |    3.90e-13 |    2.95e-14 |    5.97e-14 |    7.98e-14 |  0    0    0 
%   5.37807222799613368985795503 |     5.92e-13 |    9.08e-15 |    1.57e-12 |    8.72e-14 |    4.60e-14 |  0    0    0 
%   5.44699842002785581489333857 |     1.47e-10 |    4.42e-15 |    2.31e-14 |    2.48e-14 |    4.42e-15 |  0    0    1 
% ---------------------------------------------------------------------------------------------------------------------

% We use matrices so large that we can not use threepareigs. We use threepareigs_jd 
% which calls threepareig to solve the smaller projected three-parameter.
% We compute twice as much eigenvalues with higher accuracy and faster as 
% in demo_ellipsoidwave. In the end we use computation in multiple
% precision to refine the results.

is_numeric_type_supported('mp'); % check if MCT is installed

% we use more than quadruple precision
oldDigits = mp.Digits;
if oldDigits<100, mp.Digits(100); end

x0 = mp('1');
y0 = mp('1.5');
z0 = mp('2');

a2 = z0^2-x0^2;
b2 = z0^2-y0^2;
c = a2/b2;

multi_omega = mp([]); 
double_omega = mp([]); 

disp('Please be patient, this computation can take couple of minutes ...')
for N = [15 20 25 30 100]
    Momega = mp([]); 
    Momegad = mp([]); 
    Mlambda = mp([]); 
    Mmu = mp([]); 
    Meta = mp([]); 
    Mrst = mp([]); 

    neig = 8;
    RST = [0 0 0;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1; 1 1 1];
    for izb = 1:8
        fprintf('Computing first %d eigenmodes of type (r,s,t) = (%d,%d,%d) for N=%d:\n',neig,RST(izb,1),RST(izb,2),RST(izb,3),N)
        [omega,X1,X2,X3,xi1,xi2,xi3,lambda,mu,eta,omegad] = ellipsoid_eigs_jd_mp(x0,y0,z0,RST(izb,1),RST(izb,2),RST(izb,3),neig,N,N,N);
        Momega = [Momega; omega];
        Momegad = [Momegad; omegad];
        Mlambda = [Mlambda; lambda]; 
        Mmu = [Mmu; mu]; 
        Meta = [Meta; eta];
        Mrst = [Mrst; ones(neig,1)*RST(izb,:)];
    end

    [Momega, ord] = sort(Momega);
    Momegad = Momegad(ord); 
    Mlambda = Mlambda(ord); 
    Mmu = Mmu(ord); 
    Meta = Meta(ord); 
    Mrst = Mrst(ord,:);
    multi_omega = [multi_omega Momega(1:20)];
    double_omega = [double_omega Momegad(1:20)];
end
    
disp(' ')
disp('  "Exact" eigenmodes and errors for "exact" approximations obtained by Chebyshev collocation');
disp('------------------------------------------------------------------------------------------------------');
disp('  omega (N=100)                | error(N=30) | error(N=25) | error(N=20) | error(N=15) | rho sigma tau  ');
disp('------------------------------------------------------------------------------------------------------');
for j = 1:size(multi_omega,1)
   fprintf('%30.26f | %11.2e | %11.2e | %11.2e | %11.2e |  %d    %d    %d \n',...
       multi_omega(j,5), abs(multi_omega(j,4)-multi_omega(j,5)),abs(multi_omega(j,3)-multi_omega(j,5)),abs(multi_omega(j,2)-multi_omega(j,5)),abs(multi_omega(j,1)-multi_omega(j,5)),Mrst(j,1),Mrst(j,2),Mrst(j,3))
end
disp('------------------------------------------------------------------------------------------------------');

disp(' ')
disp('  "Exact" eigenmodes and errors for approximations obtained by Chebyshev collocation in double precision');
disp('---------------------------------------------------------------------------------------------------------------------');
disp('  omega (N=100)                | error(N=100) | error(N=30) | error(N=25) | error(N=20) | error(N=15) | rho sigma tau  ');
disp('---------------------------------------------------------------------------------------------------------------------');
for j = 1:size(multi_omega,1)
   fprintf('%30.26f | %12.2e | %11.2e | %11.2e | %11.2e | %11.2e |  %d    %d    %d \n',...
       multi_omega(j,5), abs(double_omega(j,5)-multi_omega(j,5)), abs(double_omega(j,4)-multi_omega(j,5)),abs(double_omega(j,3)-multi_omega(j,5)),abs(double_omega(j,2)-multi_omega(j,5)),abs(double_omega(j,1)-multi_omega(j,5)),Mrst(j,1),Mrst(j,2),Mrst(j,3))
end
disp('---------------------------------------------------------------------------------------------------------------------');

mp.Digits(oldDigits);