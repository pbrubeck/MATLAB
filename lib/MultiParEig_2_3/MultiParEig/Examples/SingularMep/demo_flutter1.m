%DEMO_FLUTTER1  Aerolastic flutter as a multiparameter eigenvalue problem
%
% This example computes the flutter point for the section model with 
% unsteady aerodynamics, where the Theodorsen function is approximated by 
% the Jones approximation. The example is from
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
% The particular example uses equation (2.3.21) from page 25, the results 
% in the (kappa, chi) form agree with Figure 4.15 from page 91.
%
% The flutter point is kappa = 0.8077 and chi = 0.2557

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

c1 = -0.2808i;
c2 = -0.01365;
c3 = -0.3455i;

A0 = M0 + G0;
A1 = c3*M0 + c3*G0 + G1 + 1/2*G2;
A2 = c2*M0 + c2*G0 + c3*G1 + c1*G2 + 1/2*G3;
A3 = c2*G1 + c2*G2 + c1*G3;
A4 = c2* G3;
A5 = -D0;
A6 = -c3*D0;
A7 = -c2*D0;
A8 = -K0;
A9 = -c3*K0;
A10 = -c2*K0;

Z2 = zeros(2);
I2 = eye(2);

% linearization as described on page 80
P = zeros(20); 
P([1 2],:) = [A5 A6 A7 Z2 Z2 A8 A9 A10 Z2 Z2];
for k = 1:9
    P((2*k+1):(2*k+2),(2*k-1):(2*k))=I2;
end

Q = zeros(20); 
Q([1 2],1:10) = [A0 A1 A2 A3 A4];
Q(11:12,19:20) = -I2;

R = -eye(20);
R(1:2,1:2) = Z2;
R(11:12,11:12) = Z2;

% the obtained two-parameter eigenvalue problem is singular and has 
% multiple chi, therefore we use fast = 0
opts = [];
opts.singular = 1; 
opts.fast = 0; 
opts.showrank = 1; % display steps of the staircase algorithm 
opts.rankeps = 1e-8; % the default value (1e-10) is too small for the right detection in the staircase algorithm
[chi,kappa] = twopareig(P,-Q,-R,conj(P),-conj(Q),-conj(R),opts);
kand = find(abs(imag(kappa))<1e-4); % we extract real candidates 
sol5 = [kappa(kand) chi(kand)]

tmpk = kappa(kand); tmpc = chi(kand);

% new filter is kappa>0
filter1 = find(tmpk>0);
tmpk = tmpk(filter1); tmpc = tmpc(filter1);

% new filter is chi>0
filter2 = find(tmpc>0);
tmpk = tmpk(filter2); tmpc = tmpc(filter2);

kappa_chi5 = real([tmpk tmpc]) % the solution

Upsilon_chi5 = [kappa_chi5(:,2)./kappa_chi5(:,1) kappa_chi5(:,2)]





