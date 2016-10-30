function [] = wave1D(N)
dt=6/N^2;
x=chebGrid(N);
a=zeros(5,1); a(end)=1;
u=zeros(N,2);
u(:,1)=HermitePsi(a,20*x);
u(:,2)=0;

figure(1);
h=plot(x,u(:,1));
axis manual;
ylim([-1,1]);

nframes=10000;
for i=1:nframes
    b1=spline(x,u(:,1),1-dt);
    b2=spline(x,u(:,1),-1+dt);
    u=solveRK4(u,dt);
    u(1,1)=b1;
    u(end,1)=b2;
    set(h, 'YData', real(u(:,1)));
    drawnow;
end
end

function u=solveRK4(u, dt)
% Time-stepping by Runge Kutta 4th order.
k1=dt*partialTime(u);
k2=dt*partialTime(u+k1/2);
k3=dt*partialTime(u+k2/2);
k4=dt*partialTime(u+k3);
u=u+(k1+2*k2+2*k3+k4)/6;
end

function v=partialTime(u)
v(:,1)=u(:,2);
v(:,2)=chebfftD2(u(:,1),1);
end