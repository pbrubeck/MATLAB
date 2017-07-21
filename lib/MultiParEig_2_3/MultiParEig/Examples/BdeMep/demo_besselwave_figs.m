%DEMO_BESSELWAVE_FIGS  Bessel wave system is solved as a two-parameter 
% eigenvalue problem with Chebyshev collocation and eigenmodes are plotted
%
% This example reconstructs eigenmodes in 
% Lew Yan Voon & Willatzen, Helmholtz equation in parabolic rotational 
% coordinates: application to wave problems in quantum mechanics and
% acoustics, Math. Comp. Sim. 65 (2004) 337--349.
%
% See also: DEMO_BESSELWAVE1, DEMO_BESSELWAVE2, DEMO_BESSELWAVE3, BESSELWAVE_MEP, TWOPAREIGS

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, Math. Comp. Sim. 65 (2004) 337--349.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 14.04.2015 : change colormap to use the same colors as in Matlab 2012
% Last revision: 8.9.2015

N = 60;
close all % closes all open figures

XiSq = 1;
EtaSq = 1;

% We compute many modes and then select the smallest sqrt(mu)
% We select only modes with lambda>0 (since some are close to zero and
% negative, we take lambda>-1)

BC = [1 1];  % Dirichlet & Dirichlet boundary conditions, change to [2 2] for Neumann & Neumann or [1 2] for Dirichlet & Neumann
shift = 10;  % makes sure that A1 and A2 are nonsingular for BC = [2 2] or BC = [1 2]
modes = [];
XR = [];
YR = [];
for p=0:5
    [A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = besselwave_mep(N,N,p,sqrt(XiSq),sqrt(EtaSq),BC);
    A1 = A1+shift*C1; 
    A2 = A2+shift*C2;
    [lambda,mu,XA,YA] = twopareigs(A1,B1,C1,A2,B2,C2,8);
    XRA = recover_bc(XA,G1,k1,r1);
    YRA = recover_bc(YA,G2,k2,r2);
    modes = [modes; p*ones(length(lambda),1) lambda sqrt(mu-shift)];
    XR = [XR XRA];
    YR = [YR YRA];
end

sub = find(real(modes(:,2))>-1);
XR = XR(:,sub);
YR = YR(:,sub);
posmodes = modes(sub,:);
[tilda,ord] = sort(real(posmodes(:,3)));
Table1 = posmodes(ord(1:20),:);  % first 20 modes

% mesh values for figures
phi = linspace(0,2*pi,N);
phi = phi(:);

matXA = (z1.^2)*cos(phi)';
matYA = (z1.^2)*sin(phi)';

figure('Position',[0 0 1200 800]) 

for fig = 1:9
   matZ1 = real((XR(:,ord(fig)).*YR(:,ord(fig)))*cos(Table1(fig,1)*phi)');
   matZ1 = matZ1 / max(max(abs(matZ1)));
   if max(max(matZ1))<0.99
       matZ1 = -matZ1;
   end
   subplot(3,3,fig);
   surf(matXA,matYA,matZ1)
   shading interp
   view(135,40)
   set(gca,'FontSize',12)
end
colormap jet  % to keep Matlab 2012 colors

figure('Position',[0 0 1200 800]) 

for fig=1:9
   matZ1 = real((XR(:,ord(fig)).*YR(:,ord(fig)))*cos(Table1(fig,1)*phi)');
   matZ1 = matZ1 / max(max(abs(matZ1)));
   if max(max(matZ1))<0.99
       matZ1 = -matZ1;
   end
   subplot(3,3,fig);
   surf(matXA,matYA,matZ1)
   shading interp
   axis off
   view(2)
   set(gca,'FontSize',12)
end
colormap jet  % to keep Matlab 2012 colors


