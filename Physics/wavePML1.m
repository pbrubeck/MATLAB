function [] = wavePML1( N )
% Solves the wave equation in 1D using a Perfectly Matched Layer for BCs.
% Uses Chebyshev spectral methods in the spatial dimension
% and classic Runge-Kutta for time evolution of second order PDE.
x=chebGrid(N);

% Layer
smax=1000;
width=12;
roi=width+1:N-width;
layer=[1:width, N-width+1:N];
global sigma;
sigma=zeros(N,1);
sigma(layer)=smax*((abs(x(layer))-x(width+1))/(1-x(width+1))).^3;

% Initial conditions
u0=exp(-40*(x).^2);
v=-u0;
w=[u0,v];

dt=6/N^2;
h=plot(x(roi), w(roi,1));
%ylim([-1,1]);

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