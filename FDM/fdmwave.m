function [] = fdmwave(N)
x=linspace(-1,1,N); x=x(:);
dx=2/(N-1);
global Dx;
Dx=[-1 0 1]/(2*dx);

% Initial conditions
u0=exp(-40*(x).^2);
v=-u0;
w=[u0, v];
% du/dt=dv/dx
% dv/dt=du/dx

% Perfectly Matched Layer
smax=400;
width=12;
roi=width+1:N-width;
layer=[1:width, N-width+1:N];
global sigma;
sigma=zeros(N,1);
sigma(layer)=smax*((abs(x(layer))+x(width+1))/(1+x(width+1))).^3;

dt=6/N^2;
figure(1);
h=plot(x(roi), w(roi,1));
ylim([-1,1]);

nframes=10000;
for i=1:nframes
    w=solveRK4(w,dt);
    w([1 end],:)=0;
    set(h,'YData',w(roi,1));
    drawnow;
end
end

function wt=partialTime(w)
global Dx;
global sigma;
wt=w;
ux=conv(w(:,1), Dx); % finite diferences in space
vx=conv(w(:,2), Dx);
wt(:,1)=vx(2:end-1)-sigma.*w(:,1);
wt(:,2)=ux(2:end-1)-sigma.*w(:,2);
end

function w=solveRK4(w, dt)
% Time-stepping by Runge Kutta 4th order.
k1=dt*partialTime(w);
k2=dt*partialTime(w+k1/2);
k3=dt*partialTime(w+k2/2);
k4=dt*partialTime(w+k3);
w=w+(k1+2*k2+2*k3+k4)/6;
end