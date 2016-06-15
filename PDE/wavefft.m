function [] = wavefft(N)
c=1;
L=1;
x=linspace(0,L,N);
[xx,yy]=ndgrid(x);
[ii,jj]=ndgrid(1:N);
w=c*hypot(ii*pi/L, jj*pi/L);

f=exp(-40*((xx-L/2+.4*L/2).^2+(yy-L/2).^2));
g=0*f;

a=4/N^2*dst(dst(f)')';
b=4/N^2*dst(dst(g)')'./w;

figure(1);
h=surf(xx, yy, f, 'EdgeColor', 'none');
shading interp; colormap(jet);
zlim([-1,1]); view(2);

frames=10000;
dt=6/N^2;
tf=frames*dt;
for t=0:dt:tf
    T=a.*cos(w*t)+b.*sin(w*t);
    u=dst(dst(T)')';
    set(h, 'ZData', u);
    drawnow;
end
end