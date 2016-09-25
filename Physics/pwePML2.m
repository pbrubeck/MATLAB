function [ ] = pwePML2( N )
% Solves the paraxial wave equation in 2D using a Perfectly Matched Layer for BCs.
% Uses Hermite spectral methods in two spatial dimensions
% and classic Runge-Kutta for z evolution of second order PDE.
global D k sig;
lambda=100;
k=2*pi/lambda;
[D, x]=hermD(N);
[xx,yy]=ndgrid(x);
dz=6/N^2;

% Layer
xl=3/4*x(end);
layer=abs(x)>xl;
roi=~layer;

sig=zeros(N,1);
sig(layer)=2/dz*((abs(x(layer))-xl)/(x(end)-xl)).^3;

% Initial conditions
u0=exp(2i*pi*xx-(xx.^2+yy.^2)/2);
vx=zeros(N);
vy=zeros(N);
psi=zeros(N);
w=cat(3,u0,vx,vy,psi);

figure(1);
h=surf(xx(roi,roi), yy(roi,roi), abs(w(roi,roi,1)).^2);
colormap(jet(256)); shading interp;
xlim([-xl,xl]); ylim([-xl,xl]); zlim([-1,1]); 
view(2); axis square; 

nframes=10000;
for i=1:nframes
    w=solveRK4(w,dz);
    w([1 end],:,:)=0;
    w(:,[1 end],:)=0;
    if(mod(i,5)==1)
        h.ZData=abs(w(roi,roi,1)).^2;
        drawnow;
    end
end
end

function wt=partialTime(w)
global D k sig;
u=w(:,:,1);
vx=D*u -dgmm(sig , w(:,:,2));
vy=u*D'-dgmm(sig', w(:,:,3));

wt=w;
wt(:,:,2)=vx;
wt(:,:,3)=vy;
wt(:,:,4)=-dgmm(sig,dgmm(sig',u))+dgmm(sig',vx)+dgmm(sig,vy);
wt(:,:,1)=1i/(2*k)*(D*vx+vy*D')-dgmm(sig,u)-dgmm(sig',u)+w(:,:,4);
end

function w=solveRK4(w, dz)
% z-stepping by Runge Kutta 4th order.
k1=dz*partialTime(w);
k2=dz*partialTime(w+k1/2);
k3=dz*partialTime(w+k2/2);
k4=dz*partialTime(w+k3);
w=w+(k1+2*k2+2*k3+k4)/6;
end

function A=dgmm(x,A)
A=bsxfun(@times, x, A);
end