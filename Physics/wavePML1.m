function [] = wavePML1( N )
% Solves the wave equation in 1D using a Perfectly Matched Layer for BCs.
% Uses Chebyshev spectral methods in the spatial dimension
% and classic Runge-Kutta for time evolution of second order PDE.
x=chebGrid(N);
dt=6/N^2;

% Layer
xl=0.95;
layer=abs(x)>xl;
roi=~layer;
global sigma;
sigma=zeros(N,1);
sigma(layer)=1/dt*((abs(x(layer))-xl)/(1-xl)).^3;

% Initial conditions
u0=sech(20*(x-0.5)).^2-0.7*sech(20*(x+0.5)).^2;
v0=sign(x).*u0;
w=[u0,v0];

figure(1);
h=plot(x(roi), w(roi,1));
ylim([-1,1]);

nframes=10000;
for i=1:nframes
    w=solveRK4(w,dt);
    w([1 end],:)=0;
    if(mod(i,2)==1)
        set(h, 'YData', real(w(roi,1)));
        drawnow;
    end
end
end

function wt=partialTime(w)
global sigma;
wt=w;
wt(:,1)=chebfftD(w(:,2),1)-sigma.*w(:,1);
wt(:,2)=chebfftD(w(:,1),1)-sigma.*w(:,2);
end

function w=solveRK4(w, dt)
% Time-stepping by Runge Kutta 4th order.
k1=dt*partialTime(w);
k2=dt*partialTime(w+k1/2);
k3=dt*partialTime(w+k2/2);
k4=dt*partialTime(w+k3);
w=w+(k1+2*k2+2*k3+k4)/6;
end