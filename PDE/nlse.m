function [] = nlse( N )
% NLSE NonLinear Schrodinger Equation
% 1i u_t + 1/2 u_xx +|u|^2 u = 0
nframes=1024;
t0=0;
tf=2*pi;
x0=3;

dt=1/N;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));

% Domain substitution x=tan(th), psi=sec(th).*u
[D,th]=chebD(N); th=pi/2*th; D=2/pi*D;
x=tan(th);

% Initial condition
switch(2)
case 1 % Peregrine soliton
psi=(1-4*(1+2i*t0)./(1+4*(x.^2+t0^2)))*exp(1i*t0);
case 2 % Kuznetsov-Ma breather
b=3; a=(1+sqrt(1+b^2))/4; w=2*sqrt(2*a-1);
psi=(1+(2*(2*a-1))./(sqrt(2*a)*cosh(w*x)+1));
end
u=cos(th).*psi;

% Linear propagator
H=-1/2*diag(cos(th).^4)*(eye(N)+D*D);
[V,L]=eig(H(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N);
U(2:end-1,3:end)=V;
nH=null(H);
U(:,1:2)=nH(:,1:2);
Q=U*diag(exp(-1i*dt/2*L))/U;

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
        u=u.*exp(1i*dt*(abs(sec(th).*u).^2));
        u=Q*u;
    end
    psi=sec(th).*u;
    udata(i,:)=psi;
    set(h, 'YData', abs(psi));
    drawnow;
end

figure(2);
id=abs(x)<=x0;
surf(x(id),t0:m*dt:tf,abs(udata(:,id)).^2);
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('|\Psi|^2');
end