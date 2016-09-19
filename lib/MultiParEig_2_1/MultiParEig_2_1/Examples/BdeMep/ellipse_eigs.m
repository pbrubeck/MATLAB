function [omega,G,F,eta,xi,h,indx,indy,modes,cptime] = ellipse_eigs(mode,alpha,beta,m,n)

%ELLIPSE_EIGS   eigenfrequencies of an elliptical membrane with a fixed boundary
%
% [omega,G,F,eta,xi,h,indx,indy,modes,cptime] = ELLIPSE_EIGS(mode,alpha,beta,m,n)
% returns the first m eigenfrequencies of an elliptical membrane with 
% axis a = alpha and b = beta. 
%
% mode = 1 : even case, 
% mode = 2 : odd case
%
% The method uses matrices of size n = [n1 n2] (default is [50 40]):
% A1, B1, C1, A2, B2, C2, corresponding to the Mathieu system 
% 
%   G''(s) + (a - 2q cos (2s))*G(s) = 0,  0 < s < pi/2 (angular eq.)
%   F''(t) - (a - 2q cosh(2t))*F(t) = 0,  0 < t < xi_0  (radial eq.),
%
% discretized as a two-parameter eigenvalue problem 
%
%   A1 x = lambda B1 x + mu C1 x    (angular eq.)
%   A2 y = lambda B2 x + mu C2 y    (radial eq.)
%
% using the Chebyshev collocation (CC).
%
% Optimal n1 and n2 depends on ratio a/b, required number of frequencies
% and accuracy. Good values (for the precision of 8 decimal places) are
% a = 2, b = 1, m = 100  -> n = [50 25]
% a = 2, b = 1, m = 200  -> n = [60 30]
% a = 2, b = 1, m = 300  -> n = [70 35]
% a = 2, b = 1, m = 400  -> n = [80 40]
% a = 4, b = 1, m = 100  -> n = [65 25]
% a = 4, b = 1, m = 200  -> n = [80 30]
% a = 4, b = 1, m = 300  -> n = [90 30]
% a = 8, b = 1, m = 100  -> n = [80 20]
% a = cosh(2), b = sinh(2), m = 100  -> n = [35 40]
%
% Output:
%   omega: vector of eigenfrequencies
%   G: matrix n1 x m with eigenvectors x
%   F: matrix n2 x m with eigenvectors y
%   eta: points from 0 to pi/2 where G is discretized
%   xi: points from 0 to xi_0 where F is discretized
%   h : sqrt(alpha^2-beta^2)
%   indx: index of G (number of zeros on [0,pi])
%   indy: index of F (number of zeros on [0,xi_0])
%   modes: version for each of the solutions (1: pi even, 2: 2pi even, 3: pi odd, 4: 2Pi odd
%   cptime: computational time
%
% To plot an eigenmode, use plot_mathieu_mode(G(:,j),F(:,j),eta,xi,h,modes(j))
% where j is 1,...,m. Examples:
%
% [omega,G,F,eta,xi,h,indx,indy,modes,cptime] = ellipse_eigs(1,2,1,100,[50 25]);
% plot_mathieu_mode(G(:,19),F(:,19),eta,xi,h,modes(19))
%
% [omega,G,F,eta,xi,h,indx,indy,modes,cptime] = ellipse_eigs(1,cosh(2),sinh(2),100,[35 40]);
% plot_mathieu_mode(G(:,19),F(:,19),eta,xi,h,modes(19))
%
% See also: PLOT_MATHIEU_MODE, DEMO_MATHIEU, MATHIEU_MODES, 

% Reference: C. I. Gheorghiu, M. E. Hochstenbach, B. Plestenjak, J. Rommes: 
% Spectral collocation solutions to multiparameter Mathieu’s system, 
% Appl. Math. Comput. 218 (2012) 11990-12000.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<5, n = [50 25]; end
if mode==1, premik=0; else premik=2; end

% how many eigenvalues do we compute for pi and 2*pi case each
k = floor(m/2) + 10;

n1 = n(1);
n2 = n(2);

tic

h = sqrt(alpha^2-beta^2);

% Pi solutions
% ----------------------------------------------------------------------
[A1,B1,C1,A2,B2,C2,eta,xi,G1,k1,r1,G2,k2,r2] = mathieu_mep(n1,n2,premik+1,alpha,beta);
[lambda,mu,X1,X2] = twopareigs(A1+5*B1,B1,C1,A2+5*B2,B2,C2,k); % we shift so that A1 and A2 are not singular
XR1 = recover_bc(X1,G1,k1,r1);
XR2 = recover_bc(X2,G2,k2,r2);
for j=1:k
    indx(j,1) = 2*count_sign_changes(XR1(:,j))+premik;
    indy(j,1) = count_sign_changes(XR2(:,j))+1;
end
omega1 = sqrt(mu)/h*2;
F = XR2;
G = XR1;

% 2*Pi solutions
% ----------------------------------------------------------------------
[A1,B1,C1,A2,B2,C2,eta,xi,G1,k1,r1,G2,k2,r2] = mathieu_mep(n1,n2,premik+2,alpha,beta); 
[lambda,mu,X1,X2] = twopareigs(A1,B1,C1,A2,B2,C2,k);
XR1 = recover_bc(X1,G1,k1,r1);
XR2 = recover_bc(X2,G2,k2,r2);
for j=1:k
    indx(k+j,1) = 2*count_sign_changes(XR1(:,j))+1;
    indy(k+j,1) = count_sign_changes(XR2(:,j))+1;
end
omega2 = sqrt(mu)/h*2;
F = [F XR2];
G = [G XR1];

% union of solutions
% ----------------------------------------------------------------------
modes = [(premik+1)*ones(k,1); (premik+2)*ones(k,1)];
[omega, ord] = sort([omega1; omega2]);
G = G(:,ord(1:m));
F = F(:,ord(1:m));
indx = indx(ord(1:m),1);
indy = indy(ord(1:m),:);
omega = omega(1:m);
modes = [(premik+1)*ones(k,1); (premik+2)*ones(k,1)];
modes = modes(ord(1:m),:);

cptime = toc;





            

                
                