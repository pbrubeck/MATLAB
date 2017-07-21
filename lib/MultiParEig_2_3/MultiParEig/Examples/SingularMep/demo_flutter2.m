%DEMO_FLUTTER1  Aerolastic flutter as a multiparameter eigenvalue problem
%
% This example computes the flutter point for the section model with 
% unsteady aerodynamics, where the Theodorsen function is approximated by 
% the fractional approximation. The example is from
% 
% Arion Pons: Aerolastic flutter as a multiparameter eigenvalue problem,
% MSc thesis, Department of Mechanical Engineering, University of 
% Canterbury, 2015.
% 
% The main idea is: we have a complex matrix value function A(hi,p), where
% hi and p are complex and real parameter, respectively. We want to find a 
% p such that the eigenvalue hi has zero imaginary part. We take the 
% adjungate of the equation, assume that hi and p are real, and thus obtain 
% a nonlinear two-parameter eigenvalue problem. We linearize it into a
% polynomial two-parameter eigenvalue problem and solve it. From the 
% obtained solution we extract eigenvalues (hi,p) where both hi and p are 
% real, this is the solution of the original problem. For more details, see
% the above MSc thesis.
% 
% The particular example uses equation (2.3.28) from page 26, the results 
% in the (kappa, chi) form agree with Figure 4.11 from page 86.
%
% The flutter point is kappa = 0.88895 and chi = 0.2813

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 08.11.2016

kappa = 20;             % mass ratio
r = 0.4899;             % radius of gyration
omega_h = 0.5642;       % bending natural frequency
omega_theta = 1.4105;   % torsional natural frequency
r_theta = -0.1;         % static imbalance
a = -0.2;               % centre of mass location
d_h = 1;                % bending damp. coeff.
d_omega = 1;            % torsional damp. coeff
zeta_h = 0.01*1.4105;        % bending damp. ratio
zeta_theta = 0.01*2.3508;    % torsional damp. ratio

G0 = 1/kappa*[1 -a; -a  1/8+a^2];
G1 = 1/kappa*[0 -1i; 0 -1i*(1/2-a)];
G2 = 1/kappa*[-2i  -2i*(1/2-a); 2i*(1/2+a)  2i*(1/4-a^2)];
G3 = 1/kappa*[0  -2; 0  2*(1/2+a)];

M0 = [1  -r_theta  ; -r_theta   r^2];
D0 = [2i*zeta_h*omega_h   0; 0  2i*r^2*zeta_theta*omega_theta];
K0 = [omega_h^2  0; 0  r^2*omega_theta^2];

beta = 5/6;
F = 2.19;

A0 = 2*(1i)^beta*F*(M0 + G0);
A1 = M0 + G0;
A2 = (1i)^beta*F*(2*G1+G2);
A3 = G1 + G2;
A4 = (1i)^beta*F*G3;
A5 = G3;
A6 = -2*(1i)^beta*F*D0;
A7 = -D0;
A8 = -2*(1i)^beta*F*K0;
A9 = -K0;

Z2 = zeros(2);
I2 = eye(2);

% linearization as on p. 82
P = zeros(48); 
P(1:2,13:14) = A8;
P(1:2,23:24) = A9;
P(3:48,1:46) = -eye(46);

Q = zeros(48); 
Q(1:2,1:14) = [A0 A1 A2 A3 A4 A5 A6];
Q(1:2  ,23:24) = A7;
Q(3:4  ,13:14) = I2; % k^17 xi = xi*k^17
Q(5:6  ,23:24) = I2; % k^12 xi = xi*k^12
Q(7:8  ,25:26) = I2; % k^11 xi = xi*k^11
Q(9:10 ,35:36) = I2; % k^6 xi = xi*k^6
Q(11:12,37:38) = I2; % k^5 xi = xi*k^5
Q(13:14,47:48) = I2; 

R = zeros(48);
R(15:48,15:48) = eye(34);

% the problem is singular and has eigenvalues with multiple chi
opts = [];
opts.singular = 1; 
opts.fast = 0; 
opts.showrank = 1; 
opts.rrqr = 1;  % we use faster rank revealing QR instead of more accurate SVD
[chi,kappa] = twopareig(P,-Q,-R,conj(P),-conj(Q),-conj(R),opts);
kand = find(abs(imag(kappa))<1e-8);
realsol = [kappa(kand).^6 chi(kand)]

tmpk = kappa(kand); tmpc = chi(kand);
% new filter is kappa>0
filter1 = find(tmpk>0);
tmpk = tmpk(filter1); tmpc = tmpc(filter1);

% new filter is chi>0
filter2 = find(tmpc>0);
tmpk = tmpk(filter2); tmpc = tmpc(filter2);

kappa_chi6 = real([tmpk.^6 tmpc])

upsilon_hi6 = [kappa_chi6(:,2)./kappa_chi6(:,1) kappa_chi6(:,2)]




