function [] = wavefft(N)
%WAVEFFT Summary of this function goes here
%   Detailed explanation goes here
c=1;
L=1;
x=linspace(0,L,N);
[xx,yy]=meshgrid(x);
[ii,jj]=meshgrid(1:N);
w=c*hypot(ii*pi/L, jj*pi/L);

f=zeros(N,N);
f(hypot(xx-L/2,yy-L/2)<L/4)=1;
g=0*f;


a=4/N^2*dst(dst(f)')';
b=4/N^2*dst(dst(g)')'./w;

figure(1);
h=surf(xx, yy, f, 'EdgeColor', 'none');
shading interp
colormap(jet);
zlim([-2,2]);
set(gca, 'zlimmode','manual');


frames=512;
tf=10;
dt=tf/frames;
for t=0:dt:tf
    T=a.*cos(w*t)+b.*sin(w*t);
    u=dst(dst(T)')';
    set(h, 'ZData', u);
    drawnow;
end

end

