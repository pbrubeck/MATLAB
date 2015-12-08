function [Hol] = hologram(Nx, Ny, m, k, p, l, w)
% Generates Laguerre Gauss hologram
% Nx, Ny : number of pixels in one dimensions
% m      : number of grids
% k      : grid orientation
% p, l   : Laguerre parameters (radial order, topological charge)
% w      : waist [-1,1]^2
x=linspace(-Nx/Ny,Nx/Ny,Nx);
y=linspace(-1,1,Ny);
[xx, yy]=meshgrid(x, y);
th=atan2(yy, xx);
rr2=xx.^2+yy.^2;
qq=rr2/w^2;
L=Laguerre(p+1, abs(l));
Lag=Horner(L(end,:), 2*qq(1:end));
Lag=reshape(Lag, [Ny, Nx]);
E=(2*qq).^(abs(l)/2).*exp(-qq+1i*l*th).*Lag;
Phi=angle(E); I=E.*conj(E);
k=k/norm(k);
plane=k(1)*xx+k(2)*yy;
Hol=I.*mod(Phi+m*plane, 2*pi);
a=min(Hol(:)); b=max(Hol(:));
Hol=uint8(255*(Hol-a)/(b-a));
figure(1); image(x,y,Hol);
map=gray(256); colormap(map); colorbar();
%fullscreen(cat(3,Hol,Hol,Hol), 1);
%imwrite(Hol, map, 'Test.png');
end