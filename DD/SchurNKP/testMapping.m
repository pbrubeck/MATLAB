% Sample gird
N=128;
x=linspace(-1,1,N)';
y=linspace(-1,1,N)';
[xx,yy]=ndgrid(x,y);

% Vertices (EN, WN, ES, WS)
z=[(1+1i), (-1+1i)*4, (1-1i), (-2-1i)];
% Radius of curvautre (E, W, N, S)
r=[-2,-6,-inf,4];

F = curvedquad(z,r);
[jac,G11,G12,G22] = diffgeom(F,x,y); 
zz = F(xx,yy);

figure(1);

subplot(2,2,1);
surf(real(zz), imag(zz), jac);
title('Jacobian');
colormap(jet(256));
shading interp;
axis square;
%camlight;
view(2);

subplot(2,2,2);
surf(real(zz), imag(zz), G11);
title('Metric G11');
colormap(jet(256));
shading interp;
axis square;
%camlight;
view(2);

subplot(2,2,3);
surf(real(zz), imag(zz), G12);
title('Metric G12');
colormap(jet(256));
shading interp;
axis square;
%camlight;
view(2);

subplot(2,2,4);
surf(real(zz), imag(zz), G22);
title('Metric G22');
colormap(jet(256));
shading interp;
axis square;
%camlight;
view(2);