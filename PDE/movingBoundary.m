function [] = movingBoundary( n )

b0=1;
e=0.4;
w=0.05/e;

a =@(t) (b0+e*sin(w*t))/2;
da=@(t) e*w*cos(w*t)/2;
b =@(t) (b0+e*sin(w*t))/2;
db=@(t) e*w*cos(w*t)/2;

rd=[1,n];
[D,xi]=chebD(n);

aa=1;
bb=0;
E=eye(n);
B=diag(aa)*E(rd,:)+diag(bb)*D(rd,:);
P=eye(n)-B'/(B*B')*B;
D2=P*D*D;
D=P*D;

function dw=partialT(w,t)
    u=real(w);
    v=imag(w);
    F=(da(t)*xi+db(t))*a(t);
    G=-(da(t)*xi+db(t)).^2+a(t)^2;
    dw=v+1i*(D2*u+2*F.*(D*v))./G;
end

function w=solveRK4(w,t)
    k1=dt*partialT(w,      t);
    k2=dt*partialT(w+k1/2, t+dt/2);
    k3=dt*partialT(w+k2/2, t+dt/2);
    k4=dt*partialT(w+k3,   t+dt  );
    w=w+(k1+2*k2+2*k3+k4)/6;
end

t=0;
x=a(t)*xi+b(t);
u0=exp(-100*(x-1/2).^2);
v0=zeros(size(u0));
w=u0+1i*v0;
y=real(w);

tf=100;
dt=4/n^2;
nframes=tf/dt;
figure(1);
h=plot(x, y, '.-', 'LineWidth', 1);
xlim([0,2]);
ylim([-2; 2]*max(abs(real(w))));
axis manual;
drawnow;

last=0;

for i=1:nframes
    t=t+dt;
    w=solveRK4(w,t);
    if(t-last>0.05 || i==nframes)
        y=real(w);
        x=a(t)*xi+b(t);
        set(h,'XData',x);
        set(h,'YData',y);
        drawnow;
        last=t;
    end
end

end

