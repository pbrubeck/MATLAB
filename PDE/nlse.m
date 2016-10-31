function [] = nlse( N )
% NLSE NonLinear Schrodinger Equation
% 1i u_t + 1/2 u_xx +|u|^2 u = 0
t0=-3;
tf=3;
dt=1/(N);
x0=3;

% Domain substitution x=tan(th), psi=sec(th).*u
[D,th]=chebD(N); th=pi/2*th; D=2/pi*D;
x=tan(th);

% Initial condition
psi=(1-4*(1+2i*t0)./(1+4*(x.^2+t0^2)))*exp(1i*t0);
u=cos(th).*psi;

% Linear propagator, exponential map
H=-1/2*diag(cos(th).^4)*(eye(N)+D*D);
[V,L]=eig(H(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N);
U(:,1:2)=null(H);
U(2:end-1,3:end)=V;
Q=U*diag(exp(-1i*dt/2*L))/U;

figure(1);
h=plot(x, abs(psi), 'b', 'LineWidth', 2);
xlim([-x0,x0]); ylim([0,3]); axis manual;
xlabel('x'); title('|\Psi|');
drawnow;

nframes=N;
udata=zeros(nframes,N);
udata(1,:)=abs(psi);
for i=2:nframes
    for t=0:dt:(tf-t0)/(nframes+1)
        u=Q*u;
        u=u.*exp(1i*dt*(abs(sec(th).*u).^2));
        u=Q*u;
    end
    psi=sec(th).*u;
    udata(i,:)=psi;
    set(h, 'YData', abs(psi));
    drawnow;
end

figure(2);
id=abs(x)<=3;
surf(x(id),linspace(t0,tf,nframes),abs(udata(:,id)).^2);
colormap(jet(50)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('|\Psi|^2');
end