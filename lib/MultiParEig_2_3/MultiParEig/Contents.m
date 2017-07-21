% Multiparameter Eigenproblems (MultiParEig) Toolbox 
% Version 2.3             03-Jan-2017
% 
% Toolbox for multiparameter eigenvalue problems
% 
% Toolbox contains numerical methods for multiparameter eigenvalue problems 
% 
% ------------------------------------------------------------------------
% Acknowledgements
% ------------------------------------------------------------------------
% 
% This is a joined work with Andrej Muhiè, who wrote part of the code, 
% among other things he contributed the staircase algorithm for a singular 
% multiparameter eigenvalue problem. If you use this toolbox to solve a 
% singular multiparameter eigenvalue problem, please cite the reference: 
% A. Muhiè, B. Plestenjak: On the quadratic two-parameter eigenvalue 
% problem and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542.
% 
% Other contributors:
% - Pavel Holoborodko, Advanpix LLC, http://www.advanpix.com/
%
% Methods chebdif_mp and lagdif_mp are based on chebdif and lagdif from
% Matlab toolbox DMSUITE by J.A.C. Weideman
% 
% ------------------------------------------------------------------------
% Multiparameter eigenvalue problems
% ------------------------------------------------------------------------
% 
% In a matrix two-parameter eigenvalue problem, which has the form
% 
% A1*x = lambda*B1*x + mu*C1*x, 
% A2*y = lambda*B2*y + mu*C2*y,
% 
% we are looking for an eigenvalue (lambda,mu) and nonzero eigenvectors 
% x,y such that the above system is satisfied. A two-parameter eigenvalue 
% problem is related to a pair of generalized eigenvalue problems
% 
% Delta1 z = lambda Delta0 z,
% Delta2 z = mu Delta0 z,
% 
% where Delta0, Delta1 and Delta2 are operator determinants
% 
% Delta0 = kron(C2, B1) - kron(B2, C1)
% Delta1 = kron(C2, A1) - kron(A2, C1)
% Delta2 = kron(A2, B1) - kron(B2, A1)
% 
% and z = kron(x,y). We say that the problem is nonsingular when Delta0 is
% nonsingular. The above can straightforward be generalized to three or 
% more parameters.
% 
% In many applications a partial differential equation has to be solved on 
% some domain that allows the use of the method of separation of variables. 
% In several coordinate systems separation of variables applied to the 
% Helmholtz, Laplace, or Schrödinger equation leads to a multiparameter 
% eigenvalue problem, some important cases are Mathieu’s system, Lamé’s 
% system, and a system of spheroidal wave functions. A generic 
% two-parameter boundary value eigenvalue problem has the form
% 
% p1(x1) y1''(x1) + q1(x1) y1'(x1) + r1(x2) y1(x1) 
%    = lambda s1(x1) y1(x1) + mu s2(x2) y1(x1), 
% p2(x2) y2''(x2) + q2(x2) y2'(x2) + r2(x2) y2(x2) 
%    = lambda s2(x2) y2(x2) + mu s2(x2) y2(x2),
% 
% where x1 in [a1,b1] and x2 in [a2,b2] together with the boundary 
% conditions. Such system can be discretized into a matrix two-parameter 
% eigenvalue problem, where a good method of choice is the Chebyshev 
% collocation.
% 
% Functions in the toolbox can be used to:
% 
% - compute Delta matrices for a multiparameter eigenvalue problem
% - solve a nonsingular or singular multiparameter eigenvalue problem with 
%   arbitrary number of parameters (the limit is the overall size of the 
%   corresponding Delta matrices),
% - compute few eigenvalues and eigenvectors of a two-parameter eigenvalue 
%   problem using implicitly restarted Arnoldi or Krylov-Schur method,
% - compute few eigenvalues and eigenvectors of a two- or three-parameter
%   eigenvalue problem using the Jacobi-Davidson method or the subspace
%   iteration method
% - refine an eigenpair using the tensor Rayleigh quotient iteration   
% - discretize a two- or three-parameter boundary value eigenvalue problem 
%   with the Chebyshev collocation into a  matrix two- or three-parameter 
%   eigenvalue problem,
% - solve a quadratic two-parameter eigenvalue problem.
% 
% ------------------------------------------------------------------------
% Change log from version 2.2:
%  - most of the methods now support multiprecision 
%  - many multiprecision examples using Advanpix Multiprecision Computing 
%    Toolbox in directory multiprecision
%  - twopareigs uses twopareigs_ks by default for all systems, twopareig_ks 
%    now works for complex matrices as well
%  - chebdif and lagdif from DMSUITE are replaced by chebdif_mp and
%    lagdif_mp that support multiprecision
%  - new stopping criteria in TRQI when residual stops decreasing 
%  - sylvreal is removed, sylv now supports real and complex systems
%  - minor fixes
%  - Contents.m is changes so that MultiParEig shows in the list of
%    installed packages
%  
% ------------------------------------------------------------------------
% Change log from version 2.1:
%  - option of using a faster rank revealing qr for singular problems
%  - improved Jacobi-Davidson method for two and three parameters
%  - new report for singular problems with rank choices and gaps
%  - improved twopareig for eigenvalues with multiple lambda or mu
%  - option novector for twopareig and threepareig
%  - minor fixes
%  - additional examples (demo_flutter1, demo_flutter2)
%  - updates required by the forthcoming new release of BiRoots
% 
% ------------------------------------------------------------------------
% Change log from version 1.0:
%  - new numerical methods (Jacobi-Davidson and subspace iteration) that
%    can be applied to two- or three-parameter eigenvalue problems with 
%    large matrices
%  - support for singular problems with three or more parameters
%  - tensor Rayleigh quotient iteration generalized to three or more
%    parameters
%  - functions for the construction of Delta operator determinants
%  - lapack package is not needed anymore in Matlab 2014a or later
%  - additional examples for new functions and methods
% 
% ------------------------------------------------------------------------
% Main Matlab functions in the toolbox 
% ------------------------------------------------------------------------
% 
% Two-parameter eigenvalue problems (2EP):
% 
% - twopareig: solve a 2EP (set options to solve a singular 2EP)
% - twopareigs: few eigenvalues and eigenvectors of a 2EP using implicitly 
%               restarted Arnoldi method or Krylov-Schur method
% - twopareigs_si: subspace iteration with Arnoldi expansion for a 2EP
% - twopareigs_jd: Jacobi-Davidson method for a 2EP
% - trqi: tensor Rayleigh quotient iteration for a 2EP
% - twopar_delta: Delta matrices for a 2EP
% 
% Three-parameter eigenvalue problems (3EP):
% 
% - threepareig: solve a 3EP (set options to solve a singular 3EP)
% - threepareigs: few eigenvalues and eigenvectors of a 3EP using
%                 implicitly restarted Arnoldi method
% - threepareigs_si: subspace iteration with Arnoldi expansion for a 3EP
% - threepareigs_jd: Jacobi-Davidson method for a 3EP
% - trqi_3p: tensor Rayleigh quotient iteration for a 3EP
% - threepar_delta: Delta matrices for a 3EP
% 
% Multi-parameter eigenvalue problems (MEP):
% 
% - multipareig: solve a MEP (set options to solve a singular MEP)
% - trqi_np: tensor Rayleigh quotient iteration for a MEP
% - multipar_delta: Delta matrices for a MEP
% 
% Two and three-parameter boundary differential equations:
% 
% - bde2mep: discretizes two-parameter BDE as a two-parameter matrix 
%            pencil using the Chebyshev collocation
% - bde3mep: discretizes three-parameter BDE as a three-parameter matrix 
%            pencil using the Chebyshev collocation
% 
% Quadratic two-parameter eigenvalue problem (Q2EP):
% 
% - quad_twopareig: eigenvalues and eigenvectors of a Q2EP
% - linearize_quadtwopar: linearize Q2EP as a linear two-parameter matrix 
%                         pencil 
% Other applications:
% 
% - double_eig: lambda such that A+lambda*B has a multiple eigenvalue
% 
% See directory Examples with many demos. In particular, directory BdeMep 
% contains demo functions with numerical examples from B. Plestenjak, C.I. 
% Gheorghiu, M.E. Hochstenbach: Spectral collocation for multiparameter 
% eigenvalue problems arising from separable boundary value problems, 
% J. Comp. Phys. 298 (2015) 585-601
% 
% ------------------------------------------------------------------------
% References
% ------------------------------------------------------------------------
% 
% This is a list of references containing some of the implemented algorithms. 
% If you use the toolbox in your paper, please cite an appropriate reference.
% 
% twopareig (nonsingular), twopareigs_jd:
% M. E. Hochstenbach, T. Košir, B. Plestenjak: A Jacobi-Davidson type method 
% for the two-parameter eigenvalue problem, SIAM J. Matrix Anal. Appl. 26 
% (2005) 477-497
% 
% twopareig (singular), quad_twopareig, linearize_quadtwopar:
% A. Muhiè, B. Plestenjak: On the quadratic two-parameter eigenvalue problem 
% and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542
% 
% twopareigs, twopareigs_si:
% K. Meerbergen, B. Plestenjak: A Sylvester-Arnoldi type method for the 
% generalized eigenvalue problem with two-by-two operator determinants, 
% Numer. Linear Algebra Appl. 22 (2015) 1131-1146
% 
% bde2mep, bde3mep:
% B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral collocation 
% for multiparameter eigenvalue problems arising from separable boundary 
% value problems, J. Comput. Phys. 298 (2015) 585-601
% 
% twopareigs_jd (harmonic Ritz values):
% M. E. Hochstenbach, B. Plestenjak: Harmonic Rayleigh-Ritz for the 
% multiparameter eigenvalue problem, Electron. Trans. Numer. Anal. 29 (2008) 
% 81-96.
% 
% double_eig:
% A. Muhiè, B. Plestenjak: A method for computing all values lambda such 
% that A + lambda*B has a multiple eigenvalue, Linear Algebra Appl. 440 
% (2014) 345-359
% 
% ------------------------------------------------------------------------
% 
% MultiParEig toolbox
% Bor Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt
