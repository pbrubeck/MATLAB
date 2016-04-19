function [] = fdmwave(n)
x=linspace(-1,1,n);
dx=2/(n-1);
u=0*exp(-100*x.^2/2);
ul=u;
c=0.5;
src=floor(n/2);

h=plot(x,u);
axis manual;
ylim([-1,1]);
dt=0.001;
Dxx=[1 -2 1]/dx^2;
p=c*dt/dx;
k=floor(p);
p=p-k;
for t=0:dt:10
    uxx=conv(Dxx, u);
    temp=u;
    u=2*u-ul+(c*dt)^2*uxx(2:end-1);
    ul=temp;    
    u(1)=p*ul(2+k)+(1-p)*ul(1+k);
    u(end)=p*ul(end-1-k)+(1-p)*ul(end-k);    
    u(src)=min(t/2,1)*sin(pi*t);
    set(h,'YData',u);
    drawnow;
end

end