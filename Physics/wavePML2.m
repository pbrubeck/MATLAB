function [] = wavePML2(N)
% Solves the wave equation in 2D using a Perfectly Matched Layer for BCs.
% Uses Chebyshev spectral methods in two spatial dimensions
% and classic Runge-Kutta for time evolution of second order PDE.
x=chebGrid(N);
[xx,yy]=ndgrid(x);

% Layer
smax=1000;
width=12;
roi=width+1:N-width;
layer=[1:width, N-width+1:N];
global sigma;
sigma=zeros(N,1);
sigma(layer)=smax*((abs(x(layer))-x(width+1))/(1-x(width+1))).^3;

% Initial conditions
u0=exp(-40*((xx-.4).^2+yy.^2));
vx=u0;
vy=zeros(N);
phi=zeros(N);
w=cat(3,u0,vx,vy,phi);

dt=6/N^2;
h=surf(xx(roi,roi), yy(roi,roi), w(roi,roi,1), 'EdgeColor', 'none');
colormap(jet(256)); alpha(0.85); shading interp;
view(2); zlim([-1,1]); axis square; 

nframes=10000;
for i=1:nframes
    w=solveRK4(w,dt);
    w([1 end],:,:)=0;
    w(:,[1 end],:)=0;
    if(mod(i,2)==1)
        set(h, 'ZData', w(roi,roi,1));
        drawnow;
    end
end
end

function wt=partialTime(w)
global sigma;
wt=w;
ux=chebfftD(w(:,:,1),1);
uy=chebfftD(w(:,:,1),2);
vx=chebfftD(w(:,:,2),1);
vy=chebfftD(w(:,:,3),2);
wt(:,:,1)=(vx+vy)-dgmm(sigma, w(:,:,1))-dgmm(sigma', w(:,:,1))+w(:,:,4);
wt(:,:,2)=ux-dgmm(sigma, w(:,:,2));
wt(:,:,3)=uy-dgmm(sigma', w(:,:,3));
wt(:,:,4)=(dgmm(sigma', vx)+dgmm(sigma, vy))-dgmm(sigma', dgmm(sigma, w(:,:,1)));
end

function w=solveRK4(w, dt)
% Time-stepping by Runge Kutta 4th order.
k1=dt*partialTime(w);
k2=dt*partialTime(w+k1/2);
k3=dt*partialTime(w+k2/2);
k4=dt*partialTime(w+k3);
w=w+(k1+2*k2+2*k3+k4)/6;
end

function A=dgmm(x,A)
A=bsxfun(@times, x, A);
end