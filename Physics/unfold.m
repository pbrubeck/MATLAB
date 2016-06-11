function []=unfold(n)
tau=0.7;
[xx,yy]=meshgrid(linspace(-1,1,n));
R=50*exp(-3*(xx.^2+yy.^2));
A=mod(R,2*pi);
lap=[-1 -1 -1; -1 8 -1; -1 -1 -1];

Dx=[0 0 0; -1 0 1; 0 0 0];
Dy=[0 -1 0,; 0 0 0; 0 1 0];

L=filter2(lap,A);
Lx=filter2(Dx,A);
Ly=filter2(Dy,A);
L=0*(abs(L/max(L(:))-1)<=tau)+(abs(L/min(L(:))-1)<=tau);
Lx=(Lx>0)-(Lx<0); Lx(L==0)=0;
Ly=(Ly>0)-(Ly<0); Ly(L==0)=0;

T=tril(ones(n));
Px=Lx*T';
Py=T*Ly;

figure(1);
surf(xx,yy,A-2*pi*min(Px,Py)); shading interp
colormap(jet(256));
colorbar();

figure(2);
imagesc(L);
colormap(gray(256));
end