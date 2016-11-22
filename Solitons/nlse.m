function [] = nlse(N, beta, init, method)
% NonLinear Schrodinger Equation in one dimension
% 1i u_t + 1/2 u_xx +|u|^2 u = 0

t0=0; tf=5; 
nframes=1024;
dt=1/1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));
x0=8;

% Linear propagator, half step
if nargin==3, method='cheb'; end;
switch(method)
    case 'cheb', [x,w,D,Q]=nlsecheb(N,dt/2);
    case 'herm', [x,w,D,Q]=nlseherm(N,dt/2);
    case 'fft',  [x,w,D,Q]=nlsefft(N, dt/2);
end

switch(init)
    case 'test', u=sech(x+5).*exp(1i*2*x)+sech(x-5).*exp(-1i*2*x);
end


E=w'*(abs(linprop(D,u)).^2+abs(u).^4);
display(E);

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
        u=u.*exp(1i*beta*dt*(abs(u).^2));
        u=linprop(Q,u);
    end
    udata(i,:)=u;
    set(h, 'YData', abs(u));
    drawnow;
end

E=w'*(abs(linprop(D,u)).^2+abs(u).^4);
display(E);

t=t0:m*dt:tf;

figure(2);
roi=abs(x)<=x0;
surf(x(roi),t,abs(udata(:,roi)).^2);
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('|\Psi|^2');

figure(3);
imagesc(x(roi),t,angle(udata(:,roi)));
colormap(hsv(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('Arg(\Psi)');
end

function u=linprop(Q,u)
if isfloat(Q)
    u=Q*u;
else
    u=Q(u);
end
end