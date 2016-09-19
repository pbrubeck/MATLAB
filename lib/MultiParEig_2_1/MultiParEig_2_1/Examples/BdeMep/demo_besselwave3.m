%DEMO_BESSELWAVE3  Bessel wave system is solved as a two-parameter eigenvalue problem 
% with Chebyshev collocation
%
% This example reconstructs Figures 6 (a,b,c) in  
% Lew Yan Voon & Willatzen, Helmholtz equation in parabolic rotational 
% coordinates: application to wave problems in quantum mechanics and
% acoustics, Math. Comp. Sim. 65 (2004) 337--349.
%
% See also: DEMO_BESSELWAVE1, DEMO_BESSELWAVE2, DEMO_BESSELWAVE_FIGS, BESSELWAVE_MEP, TWOPAREIGS

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 26.08.2015 : use the main routine twopareigs 
% Last revision: 8.9.2015

C = 2;
N = 50;
XiSq = linspace(0.5,2,100);
wavespeed = 343;  % wave speed in m/s^2

% -------------------------------------------------
% Figure a) : Dirichlet & Dirichlet BC, p = 0
mu = [];
for k = 1:length(XiSq);
    EtaSq(k) = max(roots([XiSq(k) XiSq(k)^2 -C]));
    [A1,B1,C1,A2,B2,C2] = besselwave_mep(N,N,0,sqrt(XiSq(k)),sqrt(EtaSq(k)),[1 1]);
    A1=A1+10*C1; % We shift to make sure that Delta2 is invertible
    A2=A2+10*C2;
    [tilda,tmpmu] = twopareigs(A1,B1,C1,A2,B2,C2,1);
    mu(k)=sqrt(tmpmu-10);
end
Frequency1 = mu*wavespeed/(2*pi);
plot(XiSq,Frequency1,'LineWidth',2)
set(gca,'FontSize',20)

% -------------------------------------------------
% Figure b) : Neumann & Neumann BC, p = 1
mu = [];
for k = 1:length(XiSq);
    [A1,B1,C1,A2,B2,C2] = besselwave_mep(N,N,1,sqrt(XiSq(k)),sqrt(EtaSq(k)),[2 2]);
    A1=A1+10*C1; % We shift to make sure that Delta2 is invertible
    A2=A2+10*C2;
    [tilda,tmpmu] = twopareigs(A1,B1,C1,A2,B2,C2,1);
    mu(k)=sqrt(tmpmu-10);
end
Frequency2 = mu*wavespeed/(2*pi);
figure
plot(XiSq,Frequency2,'LineWidth',2)
set(gca,'FontSize',20)

% -------------------------------------------------
% Figure c) : Dirichlet & Neumann BC, p = 0
mu = [];
for k = 1:length(XiSq);
    [A1,B1,C1,A2,B2,C2] = besselwave_mep(N,N,0,sqrt(XiSq(k)),sqrt(EtaSq(k)),[1 2]);
    A1=A1+10*C1; % We shift to make sure that Delta2 is invertible
    A2=A2+10*C2;
    [tilda,tmpmu] = twopareigs_ks(A1,B1,C1,A2,B2,C2,1);
    mu(k)=sqrt(tmpmu-10);
end
Frequency3 = mu*wavespeed/(2*pi);
figure
plot(XiSq,Frequency3,'LineWidth',2)
set(gca,'FontSize',20)

