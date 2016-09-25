function [] = pwePML1( N )
% Solves the paraxial wave equation in 1D using a Perfectly Matched Layer for BCs.
% Uses Hermite spectral methods in the spatial dimension
% and classic Runge-Kutta for z evolution of second order PDE.
global Dx k sig;
lambda=100;
k=2*pi/lambda;
[Dx,x]=hermD(N);
dz=6/N^2;

% Layer
xl=3/4*x(end);
layer=abs(x)>xl;
roi=~layer;

sig=zeros(N,1);
sig(layer)=2/dz*((abs(x(layer))-xl)/(x(end)-xl)).^3;

% Initial conditions
u0=exp(2i*pi*x-x.^2/2);
v0=zeros(size(x));
w=[u0,v0];


figure(1);
h=plot(x(roi), w(roi,1));
xlim([-xl,xl]);
ylim([-1,1]);

nframes=10000;
for i=1:nframes
    w=solveRK4(w,dz);
    w([1 end],:)=0;
    if(mod(i,2)==1)
        set(h, 'YData', abs(w(roi,1)).^2);
        drawnow;
    end
end
end

function wt=partialz(w)
global Dx k sig;
wt=w;
wt(:,2)=Dx*w(:,1)-sig.*w(:,2);
wt(:,1)=1i/(2*k)*Dx*wt(:,2)-sig.*w(:,1);
end

function w=solveRK4(w, dz)
% z-stepping by Runge Kutta 4th order.
k1=dz*partialz(w);
k2=dz*partialz(w+k1/2);
k3=dz*partialz(w+k2/2);
k4=dz*partialz(w+k3);
w=w+(k1+2*k2+2*k3+k4)/6;
end