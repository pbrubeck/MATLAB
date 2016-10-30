function [] = nlse( N )
% NLSE NonLinear Schrodinger Equation
% 1i u_t + 1/2 u_xx +|u|^2 u = 0
x0=20;
[D,x]=chebD(N); x=x0*x; D=D/x0;
dt=6/N^2;

% Linear propagator, exponential map
H=-1/2*D*D;
[V,L]=eig(H(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N,N);
U(:,1:2)=null(H);
U(2:end-1,3:end)=V;
Q=U*diag(exp(-1i*dt*L))/U;

% Initial conditions
t0=-1; tf=1;
u=(1-4*(1+2i*t0)./(1+4*(x.^2+t0^2)))*exp(1i*t0);

figure(1);
h=plot(x, abs(u), 'b', 'LineWidth', 1);
ylim([0,3]);
axis manual; drawnow;

nframes=1000;
for i=1:nframes
    for t=0:dt:(tf-t0)/nframes
        u=solveRK4(Q*u,dt);
    end
    set(h, 'YData', abs(u));
    drawnow;
end
end

function wt=partialt(w)
% nonlinear operator
wt=1i*(abs(w).^2).*w;
end

function w=solveRK4(w, dt)
% time-stepping by Runge Kutta 4th order.
k1=dt*partialt(w);
k2=dt*partialt(w+k1/2);
k3=dt*partialt(w+k2/2);
k4=dt*partialt(w+k3);
w=w+(k1+2*k2+2*k3+k4)/6;
end