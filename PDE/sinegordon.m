function [] = sinegordon( N )
% Sine-Gordon Equation
% u_tt - u_xx + sin(u) = 0
nframes=1024;
t0=-10;
tf=10;
x0=0.75*(tf-t0)/2;

dt=1/1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));

% Domain substitution x=tan(th)
[D,th]=chebD(N); th=pi/2*th; D=2/pi*D;
x=tan(th);

% Initial condition
[psi1,phi1]=kink(-0.75,0,x,t0);
[psi2,phi2]=kink(0.75,0,x,t0);
psi=psi1+psi2;
phi=phi1+phi2;

% Linear propagator
Dx=diag(cos(th).^2)*D;
[V,L]=eig(Dx(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N);
U(2:end-1,3:end)=V;
nullspace=null(Dx);
U(:,1:2)=nullspace(:,1:2);
Q1=real(U*diag(exp(dt/2*L))/U);
Q2=real(U*diag(sinc(1i*dt/2*L)*dt/2)/U);
Q3=real(U*diag(exp(-dt/2*L))/U);

% 2d plot
figure(1);
h=plot(x, psi, 'b', 'LineWidth', 2);
xlim([-x0,x0]); ylim([0,8]); axis manual;
xlabel('x'); title('\Psi(x)');
drawnow;

% 3d plot
figure(2);
w=exp(1i*psi);
h3=plot3(x, real(w), imag(w), 'b', 'LineWidth', 2); hold on;

xq=linspace(-x0,x0,32);
wq=interp1(x,w,xq,'spline');
q3=quiver3(xq,0*xq,0*xq,0*xq,real(wq),imag(wq),'r','LineWidth',1,'AutoScale','off'); hold off;
xlim([-x0,x0]); ylim([-2,2]); zlim([-2,2]); axis manual;
xlabel('x'); ylabel('y'); zlabel('z'); title('y+iz=e^{i\Psi(x)}'); view(-145,52);
drawnow;

% Time propagation
udata=zeros(nframes,N);
udata(1,:)=psi;
for i=2:nframes
    for j=1:m
        psi=Q1*psi+Q2*phi;
        phi=Q3*phi;
        
        phi=phi-dt*sin(psi);
        
        psi=Q1*psi+Q2*phi;
        phi=Q3*phi;
    end
    udata(i,:)=psi;
    set(h, 'YData', psi);
    
    w=exp(1i*psi);
    wq=interp1(x,w,xq,'spline');
    set(h3, 'YData', real(w));
    set(h3, 'ZData', imag(w));
    set(q3, 'VData', real(wq));
    set(q3, 'WData', imag(wq));
    drawnow;
end

% Surf plot
figure(3);
id=abs(x)<=x0;
surf(x(id),t0:m*dt:tf,udata(:,id));
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('\Psi(x,t)');
end

function [psi,phi]=kink(c,d,x,t)
g=1/sqrt(1-c^2);
psi=4*atan(exp(g*(x-c*t)+d));
phi=-2*g*(c+1).*sech(g*(x-c*t)+d);
end

function [psi,phi]=antikink(c,d,x,t)
g=-1/sqrt(1-c^2);
psi=4*atan(exp(g*(x-c*t)+d));
phi=-2*g*(c+1).*sech(g*(x-c*t)+d);
end