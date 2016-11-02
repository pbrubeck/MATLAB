function [] = sinegordon( N )
% Sine-Gordon Equation
% u_tt - u_xx + sin(u) = 0
nframes=1024;
t0=0;
tf=10;
x0=6;

dt=1/(N*20);
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));

% Domain substitution x=tan(th), psi=sec(th).*u
[D,th]=chebD(N); th=pi/2*th; D=2/pi*D;
x=tan(th);

% Initial condition
switch(1)
case 1
v=0.1; g=1/(1-v^2); d=0;
psi=4*atan(exp(g*(x-v*t0)+d));
case 2
psi=exp(-x.^2);
end
u=[cos(th).*psi; zeros(N,1)];

% Linear propagator
Dx=diag(cos(th).^2)*D;
A=[Dx, eye(N); zeros(N), -Dx]; 
Q=expm(A*dt/2);

figure(1);
h=plot(x, abs(psi), 'b', 'LineWidth', 2);
xlim([-x0,x0]); ylim([0,4]); axis manual;
xlabel('x'); title('|\Psi|');
drawnow;

udata=zeros(nframes,N);
udata(1,:)=psi;
for i=2:nframes
    for j=1:m
        u=Q*u;
        u=solveRK4(u,dt);
        u=Q*u;
    end
    psi=sec(th).*u(1:N);
    udata(i,:)=psi;
    set(h, 'YData', abs(psi));
    drawnow;
end

figure(2);
id=abs(x)<=x0;
surfl(x(id),t0:m*dt:tf,abs(udata(:,id)).^2,'light');
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('|\Psi|^2');
end

function v=partialTime(u)
N=length(u)/2;
v=[0*u(1:N); -sin(u(1:N))];
end

function u=solveRK4(u, dt)
% Time-stepping by Runge Kutta 4th order.
k1=dt*partialTime(u);
k2=dt*partialTime(u+k1/2);
k3=dt*partialTime(u+k2/2);
k4=dt*partialTime(u+k3);
u=u+(k1+2*k2+2*k3+k4)/6;
end