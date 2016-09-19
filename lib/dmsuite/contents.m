% The Matlab Differentiation Matrix Suite
%
%by J.A.C. Weideman and S.C. Reddy
%
%This is a software package consisting of seventeen {\sc Matlab}
%functions for solving differential equations on bounded, periodic,
%and infinite intervals by the spectral collocation (pseudospectral)
%method.  The package includes functions for generating differentiation
%matrices of arbitrary order corresponding to Chebyshev, Hermite,
%Laguerre, Fourier, and sinc interpolants. In addition, functions
%are included for computing derivatives via the fast Fourier transform
%for Chebyshev, Fourier, and Sinc interpolants.  Auxiliary functions
%are included for incorporating boundary conditions, performing
%interpolation using barycentric formulas, and computing roots of
%orthogonal polynomials.  In the accompanying paper it is demonstrated
%how to use the package by solving eigenvalue, boundary value, and
%initial value problems arising in the fields of special functions,
%quantum mechanics, nonlinear waves, and hydrodynamic stability.
%(The paper has been published in the ACM Transactions of Mathematical
% Software,  Vol. 26, No. 4, December 2000, Pages 465-519.) 
%
%
%   ------------------------------------------------------------
%   |  We have made every effort to test these functions.      |
%   |  However, NO GUARANTEES are made about their validity.   |
%   |  Please send a email to weideman@sun.ac.za or            |
%   |  reddy@intergate.ca if you find bugs or have comments    |
%   ------------------------------------------------------------
%
%   Notes added May 2003:
%     1) The functions poldif.m and chebdif.m can probably be
%        improved marginally by the suggestions in 
%        R. Baltensperger & M.R. Trummer,  "Spectral Differencing
%        with a Twist", to appear in SIAM J Sci Comp.
%     2) Almost all our codes will break down if the order of
%        the matrix becomes too large; in the case of lagdif.m
%        it is around 100x100 or so.   If you need to use matrices 
%        this large, it begs the  question whether you should be using 
%        this particular pseudospectral method at all.
%     
%Summary of the {\sc Matlab} functions in the suite:
%--------------------------------------------------
%
%I. Differentiation Matrices (Polynomial Based)
%
%1.    poldif.m: General differentiation matrices.
%2.    chebdif.m: Chebyshev differentiation matrices. 
%3.    herdif.m: Hermite differentiation matrices. 
%4.    lagdif.m: Laguerre differentiation matrices.
%
%II. Differentiation Matrices (Non-Polynomial)
%
%1.   fourdif.m: Fourier differentiation matrices.  
%2.   sincdif.m: Sinc differentiation matrices. 
%
%III. Boundary Conditions
%
%1.   cheb2bc.m: Chebyshev 2nd derivative matrix
%                incorporating Robin conditions.
%2.    cheb4c.m: Chebyshev 4th derivative matrix
%                incorporating clamped conditions.
%
%IV. Interpolation
%
%1.    polint.m: Barycentric polynomial interpolation on
%                arbitrary distinct nodes
%2.   chebint.m: Barycentric polynomial interpolation on
%                Chebyshev nodes. 
%3.   fourint.m: Barycentric trigonometric interpolation at
%                equidistant nodes. 
%
%V. Transform-based derivatives  
%
%1. chebdifft.m: Chebyshev derivative
%2. fourdifft.m: Fourier derivative
%3.  sincdift.m: Sinc derivative 
%
%VI. Roots of Orthogonal Polynomials
%
%1.  legroots.m: Roots of Legendre polynomials.
%2.  lagroots.m: Roots of Laguerre polynomials. 
%3.  herroots.m: Roots of Hermite polynomials.  
%
%VII. Examples
%
%1.     cerfa.m: Function file for computing the complementary 
%                error function.  Boundary condition (a) is used.
%2.     cerfb.m: Same as cerfa.m but boundary condition (b) is used.
%3.   matplot.m: Script file for plotting the characteristic curves 
%                of Mathieu's equation.
%4.       ce0.m: Function file for computing the Mathieu cosine
%                elliptic function.
%5.     sineg.m: Script file for solving the sine-Gordon equation.
%6.     sgrhs.m: Function file for computing the right-hand side of 
%                the sine-Gordon system.
%7.    schrod.m: Script file for computing the eigenvalues of the 
%                Schr\"odinger equation.
%8.    orrsom.m: Script file for computing the eigenvalues of the 
%                Orr-Sommerfeld equation.
