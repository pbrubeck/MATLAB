function [] = NavierStokes(N)
N([1 2])=N;

global Dx Dxx Dy Dyy nu F;
[Dx,x]=chebD(N(1)); Dxx=Dx*Dx;
[Dy,y]=chebD(N(2)); Dyy=Dy*Dy;
[xx,yy]=ndgrid(x,y);
[xq,yq]=meshgrid(linspace(-1,1,32));


nu=0.1;
zz=xx+1i*yy;
F=-1000*(zz.^3)./abs(zz.^3);
F=F-solenoid(F);
rho=exp(-10*abs(zz-0.1).^2);

u=-250i*(xx+2i*yy).*exp(-20*abs(zz).^2);
u=solenoid(u);


figure(1); clf; hold on;
uq=interp2(xx',yy',u.',xq,yq); uq=uq./abs(uq);
hq=quiver(xq, yq, real(uq), imag(uq), 'w');
hs=surf(xx,yy,-1+0*rho,rho); shading interp; colormap(jet(256));
axis equal; xlim([-1.1 1.1]); ylim([-1.1 1.1]);
hold off;

dt=0.1/prod(N);
nframes=10000;
for i=1:nframes
    % Velocity update
    u=velocityRK4(u, dt);
    u=solenoid(u);
    u([1 end],:)=0;
    u(:,[1 end])=0;
    
    % Density update
    rho=densityRK4(rho, u, dt);
    rho(rho<0)=0;
    rho([1 end],:)=0;
    rho(:,[1 end])=0;
    
    % Plot
    uq=interp2(xx',yy',u.',xq,yq); uq=uq./abs(uq);
    hq.UData=real(uq);
    hq.VData=imag(uq);
    hs.CData=rho;
    drawnow;
end

end

function [u]=solenoid(u)
global Dx Dxx Dy Dyy;
divu=Dx*real(u)+imag(u)*Dy';
phi=zeros(size(u));
phi(2:end-1,2:end-1)=sylvester(Dxx(2:end-1,2:end-1), Dyy(2:end-1,2:end-1)', -divu(2:end-1,2:end-1));
u=u+(Dx*phi)+1i*(phi*Dy');
end

function [a]=accel(u)
global Dx Dxx Dy Dyy nu F;
a=-real(u).*(Dx*u)-imag(u).*(u*Dy')+nu*(Dxx*u+u*Dyy')+F;
end

function [rho_dot]=flow(rho, u)
global Dx Dxx Dy Dyy nu;
rho_dot=-real(u).*(Dx*rho)-imag(u).*(rho*Dy')+nu*(Dxx*rho+rho*Dyy');
end

function [u]=velocityRK4(u, dt)
k1=dt*accel(u     );
k2=dt*accel(u+k1/2);
k3=dt*accel(u+k2/2);
k4=dt*accel(u+k3  );
u=u+(k1+2*k2+2*k3+k4)/6;
end

function [rho]=densityRK4(rho, u, dt)
k1=dt*flow(rho     , u);
k2=dt*flow(rho+k1/2, u);
k3=dt*flow(rho+k2/2, u);
k4=dt*flow(rho+k3  , u);
rho=rho+(k1+2*k2+2*k3+k4)/6;
end
