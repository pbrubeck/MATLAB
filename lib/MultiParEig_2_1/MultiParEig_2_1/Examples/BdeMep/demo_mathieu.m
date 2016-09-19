%DEMO_MATHIEU  Solves Mathieu's system as a two-parameter eigenvalue problem
%and plot two eigenmodes
%
% This example plots Figures 5 and 6 in: C. I. Gheorghiu, M. E. Hochstenbach, 
% B. Plestenjak, J. Rommes: Spectral collocation solutions to multiparameter 
% Mathieu’s system, Appl. Math. Comput. 218 (2012) 11990-1200.
% 
% See also: ELLIPSE_EIGS, PLOT_MATHIEU_MODE, MATHIEU_MODES

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

% we compute first 300 eigenmodes for ellipse with a=4 and b=1
a = 4;
b = 1;
[omega,G,F,eta,xi,h,indx,indy,modes,cptim] = ellipse_eigs(1,a,b,300,[120 40]); 

% we plot eigenmode 298 (Figure 5)
ind = 298;
plot_mathieu_mode(G(:,ind),F(:,ind),eta,xi,h,modes(ind));
title(sprintf('Even eigenmode no. %d of type (%d,%d) for (a,b)=(%d,%d) with omega = %12.8f',...
    ind,indx(ind),indy(ind),a,b,omega(ind)))

figure
% we plot eigenmode 300 (Figure 6)
ind = 300;
plot_mathieu_mode(G(:,ind),F(:,ind),eta,xi,h,modes(ind));
title(sprintf('Even eigenmode no. %d of type (%d,%d) for (a,b)=(%d,%d) with omega = %12.8f',...
    ind,indx(ind),indy(ind),a,b,omega(ind)))

