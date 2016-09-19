function plot_mathieu_mode(g,f,eta,xi,h,ver)

%PLOT_MATHIEU_MODE  plots an eigenmode of an elliptical membrane, which was
% obtained by the method ellipse_eigs. See ellipse_eigs for additional information.
%
% See also: ELLIPSE_EIGS, DEMO_MATHIEU

% Reference: C. I. Gheorghiu, M. E. Hochstenbach, B. Plestenjak, J. Rommes: 
% Spectral collocation solutions to multiparameter Mathieu’s system, 
% Appl. Math. Comput. 218 (2012) 11990-12000.

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

% Angular part
tmpx = xi(end:-1:1); % -> [0,beta_0]
vecX = tmpx;
tmp = f; 
f2 = tmp(end:-1:1);

tmpy = eta(end:-1:1); % -> [0,pi/2]
vecY = [tmpy; tmpy+pi/2; tmpy+pi; tmpy+3*pi/2];
tmp = g;

% Radial part
if ver == 1 % pi even
    f1 = [tmp(end:-1:1); tmp; tmp(end:-1:1); tmp];
elseif ver == 2 % 2*pi even
    f1 = [tmp(end:-1:1); -tmp; -tmp(end:-1:1); tmp];
elseif ver == 3 % pi odd
    f1 = [tmp(end:-1:1); -tmp; tmp(end:-1:1); -tmp];
elseif ver == 4 % 2*pi odd
    f1 = [tmp(end:-1:1); tmp; -tmp(end:-1:1); -tmp];
end

matX = h*cosh(vecX)*cos(vecY)';
matY = h*sinh(vecX)*sin(vecY)';
rx = max(max(abs(matX)));
ry = max(max(abs(matY)));
matZ = f2*f1';
matZ = min(rx,ry)*matZ/max(max(abs(matZ)));
surf(matX,matY,0.5*matZ)
shading interp
axis equal

