function [] = sinegordon( N )
% Sine-Gordon Equation
% u_tt - u_xx + sin(u) = 0
nframes=1024;
t0=-6;
tf=18;
x0=6;

dt=1/(N);
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));

% Domain substitution x=tan(th)
[D,th]=chebD(N); th=pi/2*th; D=2/pi*D;
x=tan(th);

% Initial condition
switch(1)
case 1
c=-0.4; g=1/(1-c^2); d=0;
psi=4*atan(exp(g*(x-c*t0)+d));
phi=-2*g*(c+1).*sech(g*(x-c*t0)+d);
end

% Linear propagator
Dx=diag(cos(th).^2)*D;
[V,L]=eig(Dx(2:end-1,2:end-1),'vector');
L=[0;0;L]; 
U=zeros(N);
U(2:end-1,3:end)=V;
nullspace=null(Dx);
U(:,1:2)=nullspace(:,1:2);
Q1=U*diag(exp(dt/2*L))/U;
Q2=U*diag(sinc(1i*dt/2*L)*dt/2)/U;
Q3=U*diag(exp(-dt/2*L))/U;

Q1=real(Q1);
Q2=real(Q2);
Q3=real(Q3);

figure(1);
h=plot(x, abs(psi), 'b', 'LineWidth', 2);
xlim([-x0,x0]); ylim([0,8]); axis manual;
xlabel('x'); title('|\Psi|');
drawnow;

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
    set(h, 'YData', abs(psi));
    drawnow;
end

figure(2);
id=abs(x)<=x0;
surf(x(id),t0:m*dt:tf,abs(udata(:,id)));
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('|\Psi|');
end