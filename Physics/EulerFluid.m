function [] = EulerFluid(N)
phi=2*pi*(1:N)/N;
[Dr,r]=operators(N);
xx=r*cos(phi);
yy=r*sin(phi);

[ff,rr]=meshgrid(phi,r);
c=diag(cos(phi));
s=diag(sin(phi));

u=zeros(N,N,2);
[ux,uy]=pol2Cart(u(:,:,1),u(:,:,2),c,s);

dt=2/N^2;
figure(1);
h=quiver(xx, yy, ux, uy);
nframes=2000;
for i=1:nframes
    u=solveRK4(u,Dr,r,dt);
    u(1,:,:)=0;
    if(mod(i,2)==1)
        [ux,uy]=pol2Cart(u(:,:,1),u(:,:,2),c,s);
        set(h, 'UData', ux);
        set(h, 'VData', uy);
        drawnow;
    end
end
end

function [ux,uy]=pol2Cart(ur, uf, c, s)
ux=ur*c-uf*s;
uy=ur*s+uf*c;
end

function [Dr, r]=operators(N)
[Dr, r]=chebD(2*N);
Dr=Dr(1:N,1:N)+Dr(1:N,end:-1:N+1);
r=r(1:N);
end

function [m]=materialD(u, Dr, r)
ur=u(:,:,1);
uf=u(:,:,2);
m=zeros(size(u));
m(:,:,1)=ur.*(Dr*ur)+diag(1./r)*(uf.*(fftD(ur,2,1)-uf));
m(:,:,2)=ur.*(Dr*uf)+diag(1./r)*(uf.*(fftD(uf,2,1)+ur));
end

function [v]=timeD(u, Dr, r)
N=size(u,2);
F=zeros(size(u));
F(:,:,1)=-0.2*repmat(1-r,[1,size(u,2)]);
F(:,:,2)=repmat(1-r,[1,size(u,2)]);
v=-materialD(u, Dr, r)+F;
end

function [u]=solveRK4(u, Dr, r, dt)
k1=dt*timeD(u,      Dr, r);
k2=dt*timeD(u+k1/2, Dr, r);
k3=dt*timeD(u+k2/2, Dr, r);
k4=dt*timeD(u+k3,   Dr, r);
u=u+(k1+2*k2+2*k3+k4)/6;
end