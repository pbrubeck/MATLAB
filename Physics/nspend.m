function [] = nspend(N,g,k,nsteps)
T=20*pi/sqrt(g);
x0=[zeros(N,1);pi/2*ones(N,1)];
v0=zeros(2*N,1);
u=[x0,v0];
figure(1);
r=u(1:N,1);
th=u(N+1:2*N,1);
z=[0;-1i*cumsum((1+r).*exp(1i*th))/N];
hp=plot(real(z),imag(z),'.-b');
xlim([-2,2]);
ylim([-2,1]);
axis manual;
pbaspect([1 1 1]);
t=linspace(0,T,nsteps);
h=t(2)-t(1);
for i=2:nsteps
    u=rk4step(h,g,k,u);
    r=u(1:N,1);
    th=u(N+1:2*N,1);
    z=[0;-1i*cumsum((1+r).*exp(1i*th))/N];
    set(hp,'XData',real(z));
    set(hp,'YData',imag(z));
    drawnow;
end
end
function [u]=rk4step(h,g,k,u)
k1=h*partialT(g,k,u);
k2=h*partialT(g,k,u+k1/2);
k3=h*partialT(g,k,u+k2/2);
k4=h*partialT(g,k,u+k3);
u=u+(k1+2*k2+2*k3+k4)/6;
end
function [ut]=partialT(g,k,u)
N=size(u,1)/2;
r=u(1:N,1);
vr=u(1:N,2);
th=u(N+1:2*N,1);
vt=u(N+1:2*N,2);
[ii,jj]=ndgrid(1:N);
[ti,tj]=ndgrid(th);
C=cos(tj-ti).*(N+min(1-ii,1-jj));
S=sin(tj-ti).*(N+min(1-ii,1-jj));
D=diag([ones(N,1); 1+r]);
A=D*[C,-S;S,C]*D;
B=D*[S,C;-C,S]*D;
f=[g*((N:-1:1)'.*cos(th))-k*r; -g*((N:-1:1)'.*sin(th).*(1+r))];
w=[2*vr.*vt; vt.^2];
ut=[u(:,2), A\(f+B*w)];
end