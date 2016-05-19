function [] = fdmwave(n)
x=linspace(-1,1,n); x=x(:);
dx=2/(n-1);
a=zeros(5,1); a(end)=1;
u=HermitePsi(a,20*x);
ul=u;
c=0.5;
src=floor(n/2);

h=plot(x,u);
axis manual;
ylim([-1,1]);
dt=0.001;
Dxx=[1 -2 1]/dx^2;
for t=0:dt:10
    b1=spline(x,u,1-c*dt);
    b2=spline(x,u,-1+c*dt);
    uxx=conv(Dxx, u);
    temp=u;
    u=2*u-ul+(c*dt)^2*uxx(2:end-1);
    ul=temp;    
    u(1)=b1;
    u(end)=b2;    
    %u(src)=0.5*sin(pi*t);
    set(h,'YData',u);
    drawnow;
end

end