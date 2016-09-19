function ener = energy_qd(xi2)

%ENERGY_QD ground state energy
%
% ener = ENERGY_QD(xi2) is an auxiliary function that returns ground state 
% energy from Figure 5 in Lew Yan Voon & Willatzen, Helmholtz equation in 
% parabolic rotational coordinates: application to wave problems in quantum 
% mechanics and acoustics, Math. Comp. Sim. 65 (2004) 337--349.
%
% SEE DEMO_BESSELWAVE2, BESSELWAVE_MEP, TWOPAREIGS

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 11.09.2014

V = pi/2*70^3;
C = 4*V/pi;
N = 50;

eta2 = max(roots([xi2 xi2^2 -C]));
    
[A1,B1,C1,A2,B2,C2] = besselwave_mep(N,N,0,sqrt(xi2)/10,sqrt(eta2)/10,[1 1]);
[lambda2,mu2,X1,Y1] = twopareigs(A1,B1,C1,A2,B2,C2,1);

hbar    = 6.58211928e-16;          % in eV*s
mstar   = 0.067*0.510998910*1e6;   % eV/c^2 
clight  = 299792458;               % speed of light c in m/s^2
scaling = 1e8; % 1m / 100Angstrom
mili    = 1e3; % results are in meV

ener = mu2*hbar^2*clight^2/(2*mstar)*mili*scaling^2;
