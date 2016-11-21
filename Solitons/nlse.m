function [] = nlse(N, init, method)
% NonLinear Schrodinger Equation in one dimension
% 1i u_t + 1/2 u_xx +|u|^2 u = 0

t0=-pi; tf=pi; 
nframes=1024;
dt=1/1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));
x0=3;

% Linear propagator, half step
if nargin==2, method='cheb'; end;
switch(method)
    case 'cheb', [Q,x]=nlsecheb(N,dt/2);
    case 'herm', [Q,x]=nlseherm(N,dt/2);
    case 'fft',  [Q,x]=nlsefft(N, dt/2);
end
if(init==0)
    u=peregrine(x,t0);
else
    u=kmbreather(init,x,t0);
end

figure(1);
h=plot(x, abs(u), 'b', 'LineWidth', 2);
xlim([-x0,x0]); ylim([0,4]); axis manual;
xlabel('x'); title('|\Psi|');
drawnow;

udata=zeros(nframes,length(u));
udata(1,:)=u;
for i=2:nframes
    for j=1:m
        u=linprop(Q,u);
        u=u.*exp(1i*dt*(abs(u).^2));
        u=linprop(Q,u);
    end
    udata(i,:)=u;
    set(h, 'YData', abs(u));
    drawnow;
end

figure(2);
roi=abs(x)<=x0;
surf(x(roi),t0:m*dt:tf,abs(udata(:,roi)).^2);
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('|\Psi|^2');
end

function u=linprop(Q,u)
if isfloat(Q)
    u=Q*u;
else
    u=Q(u);
end
end