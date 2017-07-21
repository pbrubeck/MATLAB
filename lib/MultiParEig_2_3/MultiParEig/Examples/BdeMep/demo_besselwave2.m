%demo_besselwave2  Bessel wave system is solved as a two-parameter eigenvalue problem 
% with Chebyshev collocation
%
% This example reconstructs Figure 5 in  
% Lew Yan Voon & Willatzen, Helmholtz equation in parabolic rotational 
% coordinates: application to wave problems in quantum mechanics and
% acoustics, Math. Comp. Sim. 65 (2004) 337--349.
%
% See also: DEMO_BESSELWAVE1, DEMO_BESSELWAVE3, DEMO_BESSELWAVE_FIGS, BESSELWAVE_MEP, TWOPAREIGS

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 26.08.2015 : use the main routine twopareigs 
% Last revision: 8.9.2015

V = pi/2*70^3;
C = 4*V/pi;
N = 60;

tic

XiSq = linspace(8,240,100);
mu = [];
for k = 1:length(XiSq);
    EtaSq(k) = max(roots([XiSq(k) XiSq(k)^2 -C]));
    [A1,B1,C1,A2,B2,C2] = besselwave_mep(N,N,0,sqrt(XiSq(k))/10,sqrt(EtaSq(k))/10,[1 1]);
    [tilda,tmpmu] = twopareigs(A1,B1,C1,A2,B2,C2,1);
    mu(k) = tmpmu;
end

% physical constants and scaling
hbar    = 6.58211928e-16;          % in eV*s
mstar   = 0.067*0.510998910*1e6;   % eV/c^2 
clight  = 299792458;               % speed of light c in m/s^2
scaling = 1e8; % 1m / 100Angstrom
mili    = 1e3; % results are in meV

Energy = mu*hbar^2*clight^2/(2*mstar)*mili*scaling^2;

plot(XiSq,Energy,'LineWidth',2)
set(gca,'FontSize',20)

toc

% Find exact position of local minimums

format long e
mina = fminsearch(@energy_qd,10,optimset('TolX',1e-10))
valuea = energy_qd(mina)
minb = fminsearch(@energy_qd,150,optimset('TolX',1e-10))
valueb = energy_qd(minb)
format short e