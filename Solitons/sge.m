function [] = sge(N, t1, t2, L, init, method)
% Sine-Gordon Equation
% u_tt - u_xx + sin(u) = 0
%
% N: number of collocation points
% t1: initial time
% t2: final time
% init: kink, breather
% method: cheb, fft, herm

nframes=1024;
dt=1/1024;
m=ceil((t2-t1)/(dt*(nframes-1)));
dt=(t2-t1)/(m*(nframes-1));
Amax=4;

c=sqrt(1/2);

% Linear propagator, half step
if nargin<6, method='cheb'; end;
switch(method)
    case 'cheb', [Q,x]=sgecheb(N, dt/2);
    case 'herm', [Q,x]=sgeherm(N, dt/2);
    case 'fft',  [Q,x]=sgefft(N, dt/2);
end

% Initial condition
switch(init)
case {1, 'kink'}
    [u1,v1]=kink(c,1i*pi,x,t1);
    [u2,v2]=kink(-c,0,x,t1);
    u=real(u1+u2);
    v=real(v1+v2);
case {2, 'breather'}
    [u,v]=sgbreather(c,x,t1);
end

% 2d plot
f1=figure(1);
h=plot(x, u, 'b', 'LineWidth', 2);
xlim([-L,L]); ylim([-Amax,Amax]); axis manual;
set(gca, 'XTick', [-L,0,L]);
set(gca, 'YTick', [-Amax,0,Amax]);
xlabel('x'); title('|\Psi|');
set(gca,'FontSize',36);
drawnow;

w=-1i*exp(1i*u);
xq=linspace(-L,L,32);
wq=interp1(x,w,xq,'spline');

% 3d plot
figure(2);
hp3=plot3(x, real(w), imag(w), 'b', 'LineWidth', 2); hold on;
hq3=quiver3(xq,0*xq,0*xq,0*xq,real(wq),imag(wq),'r','LineWidth',1,'AutoScale','off'); hold off;
xlim([-L,L]); ylim([-2,2]); zlim([-2,2]); axis manual;
xlabel('x'); ylabel('y'); zlabel('z'); title('y+iz=e^{i\Psi(x)}'); view(25,20);
drawnow;

% Time propagation
inc=0; tinc=0; tframe=(t2-t1)/8;
udata=zeros(nframes,N);
udata(1,:)=u;
for i=2:nframes
    for j=1:m
        [u,v]=Q(u,v);
        v=v-dt*sin(u);
        [u,v]=Q(u,v);
        
        tinc=tinc+dt;
        if(tinc>=tframe)
            print(f1,sprintf('%s\\data\\%s%02d',pwd,init,inc),'-depsc');
            tinc=0; inc=inc+1;
        end
    end
    print(f1,sprintf('%s\\data\\%s%02d',pwd,init,inc),'-depsc');
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
print(sprintf('%s\\data\\%s%02d',pwd,init,inc),'-depsc');

% Surf plot
figure(3);
id=abs(x)<=L;
surf(x(id),t1:m*dt:t2,udata(:,id));
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('\Psi(x,t)');
end