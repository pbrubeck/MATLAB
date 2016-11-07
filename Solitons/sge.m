function [] = sge(N, method)
% Sine-Gordon Equation
% u_tt - u_xx + sin(u) = 0

t0=-10; tf=10;
nframes=1024;
dt=1/1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));

c=sqrt(1/2);
x0=6;

% Linear propagator, half step
if nargin==1, method='cheb'; end;
switch(method)
    case 'cheb', [Q,x]=sgecheb(N, dt/2);
    case 'herm', [Q,x]=sgeherm(N, dt/2);
    case 'fft',  [Q,x]=sgefft(N, dt/2);
end

% Initial condition
switch(2)
case 1
    [u1,v1]=kink(c,1i*pi,x,t0);
    [u2,v2]=kink(-c,0,x,t0);
    u=real(u1+u2);
    v=real(v1+v2);
case 2, [u,v]=sgbreather(c,x,t0);
end

% 2d plot
figure(1);
h=plot(x, u, 'b', 'LineWidth', 2);
xlim([-x0,x0]); ylim([-4*pi,4*pi]); axis manual;
xlabel('x'); title('\Psi(x)');
drawnow;

w=-1i*exp(1i*u);
xq=linspace(-x0,x0,32);
wq=interp1(x,w,xq,'spline');

% 3d plot
figure(2);
hp3=plot3(x, real(w), imag(w), 'b', 'LineWidth', 2); hold on;
hq3=quiver3(xq,0*xq,0*xq,0*xq,real(wq),imag(wq),'r','LineWidth',1,'AutoScale','off'); hold off;

xlim([-x0,x0]); ylim([-2,2]); zlim([-2,2]); axis manual;
xlabel('x'); ylabel('y'); zlabel('z'); title('y+iz=e^{i\Psi(x)}'); view(25,20);
drawnow;

% Time propagation
udata=zeros(nframes,N);
udata(1,:)=u;
for i=2:nframes
    for j=1:m
        [u,v]=Q(u,v);
        v=v-dt*sin(u);
        [u,v]=Q(u,v);
    end
    udata(i,:)=u;
    set(h, 'YData', u);
    
    w=-1i*exp(1i*u);
    wq=interp1(x,w,xq,'spline');
    set(hp3, 'YData', real(w));
    set(hp3, 'ZData', imag(w));
    set(hq3, 'VData', real(wq));
    set(hq3, 'WData', imag(wq));
    drawnow;
end

% Surf plot
figure(3);
id=abs(x)<=x0;
surf(x(id),t0:m*dt:tf,udata(:,id));
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('\Psi(x,t)');
end