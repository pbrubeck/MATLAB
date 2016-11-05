function [] = sge(N, method)
% Sine-Gordon Equation
% u_tt - u_xx + sin(u) = 0

t0=-10; tf=10;
nframes=1024;
dt=1/1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));

c=0.75;
x0=c*(tf-t0)/2;

% Linear propagator, half step
if nargin==1, method='cheb'; end;
switch(method)
    case 'cheb', [Q,x]=sgecheb(N, dt/2);
    case 'herm', [Q,x]=sgeherm(N, dt/2);
    case 'fft',  [Q,x]=sgefft(N, dt/2);
end

% Initial condition
[psi1,phi1]=kink( c,0,x,t0);
[psi2,phi2]=kink(-c,0,x,t0);
psi=psi1+psi2;
phi=phi1+phi2;

% 2d plot
figure(1);
h=plot(x, psi, 'b', 'LineWidth', 2);
xlim([-x0,x0]); ylim([0,4*pi]); axis manual;
xlabel('x'); title('\Psi(x)');
drawnow;

% 3d plot
figure(2);
w=-1i*exp(1i*psi);
h3=plot3(x, real(w), imag(w), 'b', 'LineWidth', 2); hold on;

xq=linspace(-x0,x0,32);
wq=interp1(x,w,xq,'spline');
q3=quiver3(xq,0*xq,0*xq,0*xq,real(wq),imag(wq),'r','LineWidth',1,'AutoScale','off'); hold off;
xlim([-x0,x0]); ylim([-2,2]); zlim([-2,2]); axis manual;
xlabel('x'); ylabel('y'); zlabel('z'); title('y+iz=e^{i\Psi(x)}'); view(0,90);
drawnow;

% Time propagation
udata=zeros(nframes,N);
udata(1,:)=psi;
for i=2:nframes
    for j=1:m
        [psi,phi]=Q(psi,phi);
        phi=phi-dt*sin(psi);
        [psi,phi]=Q(psi,phi);
    end
    udata(i,:)=psi;
    set(h, 'YData', psi);
    
    w=-1i*exp(1i*psi);
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