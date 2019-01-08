function [] = npend(N,g,nsteps)
T=20*pi/sqrt(g);
x0=pi/2*ones(N,1);
v0=zeros(N,1);
u=[x0,v0];
figure(1);
z=[0;cumsum(-1i*exp(1i*u(:,1)))/N];
hp=plot(real(z),imag(z),'.-b');
xlim([-2,2]);
ylim([-2,1]);
axis manual;
pbaspect([1 1 1]);
t=linspace(0,T,nsteps);
h=t(2)-t(1);
for i=2:nsteps
    u=rk4step(h,g,u);
    z=[0;cumsum(-1i*exp(1i*u(:,1)))/N];
    set(hp,'XData',real(z));
    set(hp,'YData',imag(z));
    drawnow;
%     x=cumsum(-1i*exp(1i*u(:,1)));
%     v=cumsum(exp(1i*u(:,1)).*u(:,2));
%     T=(1/2)*(v'*v);
%     U=sum(g*imag(x));
%     E=T+U;
%     display(E);
end
end
function [u]=rk4step(h,g,u)
k1=h*partialT(g,u);
k2=h*partialT(g,u+k1/2);
k3=h*partialT(g,u+k2/2);
k4=h*partialT(g,u+k3);
u=u+(k1+2*k2+2*k3+k4)/6;
end
function [ut]=partialT(g,u)
N=size(u,1);
[ii,jj]=ndgrid(1:N);
[ti,tj]=ndgrid(u(:,1));
A=cos(tj-ti).*(N+min(1-ii,1-jj));
B=sin(tj-ti).*(N+min(1-ii,1-jj));
f=-g*((N:-1:1)'.*sin(u(:,1)));
ut=[u(:,2), A\(f+B*(u(:,2).^2))];
end