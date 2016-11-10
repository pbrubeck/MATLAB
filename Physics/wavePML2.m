function [] = wavePML2(N)
% Solves the wave equation in 2D using a Perfectly Matched Layer for BCs.
% Uses Chebyshev spectral methods in two spatial dimensions
% and classic Runge-Kutta for time evolution of second order PDE.
x=chebGrid(N);
[xx,yy]=ndgrid(x);
dt=6/N^2;

% Layer
xl=0.95;
layer=abs(x)>xl;
roi=~layer;
sig=zeros(N,1);
sig(layer)=1/dt*((abs(x(layer))-xl)/(1-xl)).^3;

% Initial conditions
u0=exp(-40*((xx-.4).^2+yy.^2));
v1=u0;
v2=zeros(N);
phi=zeros(N);
w=cat(3,u0,v1,v2,phi);

function wt=partialTime(w)
wt=w;
ux=chebfftD(w(:,:,1),1);
uy=chebfftD(w(:,:,1),2);
vx=chebfftD(w(:,:,2),1);
vy=chebfftD(w(:,:,3),2);
wt(:,:,1)=(vx+vy)-dgmm(sig, w(:,:,1))-dgmm(sig', w(:,:,1))+w(:,:,4);
wt(:,:,2)=ux-dgmm(sig, w(:,:,2));
wt(:,:,3)=uy-dgmm(sig', w(:,:,3));
wt(:,:,4)=dgmm(sig', vx)+dgmm(sig, vy)-dgmm(sig', dgmm(sig, w(:,:,1)));
end

function w=solveRK4(w, dt)
k1=dt*partialTime(w);
k2=dt*partialTime(w+k1/2);
k3=dt*partialTime(w+k2/2);
k4=dt*partialTime(w+k3);
w=w+(k1+2*k2+2*k3+k4)/6;
end

figure(1);
h=surf(xx(roi,roi), yy(roi,roi), w(roi,roi,1), 'EdgeColor', 'none');
colormap(jet(256)); shading interp;
view(2); zlim([-1,1]); axis square; 

nframes=10000;
for i=1:nframes
    w=solveRK4(w,dt);
    w([1 end],:,:)=0;
    w(:,[1 end],:)=0;
    if(mod(i,5)==1)
        h.ZData=w(roi,roi,1);
        drawnow;
    end
end
end

function A=dgmm(x,A)
A=bsxfun(@times, x, A);
end