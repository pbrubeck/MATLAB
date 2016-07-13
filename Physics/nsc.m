function [] = nsc(N)
% Compressible 2D Navier-Stokes equations in strong conservative form
N([1 2])=N;

global gamma mu vel; 
gamma=2; mu=0.05;

global Dx Dxx Dy Dyy;
[Dx,x]=chebD(N(1)); Dxx=Dx*Dx;
[Dy,y]=chebD(N(2)); Dyy=Dy*Dy;
[xx,yy]=ndgrid(x,y);
[xq,yq]=meshgrid(linspace(-1,1,32));

zz=xx+1i*yy;
rho=1+0*xx;
vel=200*(zz.^4+0.4^2).*exp(-10*abs(zz).^2);
vel=solenoid(vel);
E=10+0*xx;

% Conserved state vector
w=zeros([N, 4]);
w(:,:,1)=rho;
w(:,:,2)=rho.*real(vel);
w(:,:,3)=rho.*imag(vel);
w(:,:,4)=rho.*E;
w0=w;

figure(1);
u=w(:,:,2)+1i*w(:,:,3);
uq=interp2(xx',yy',u.',xq,yq); uq=uq./abs(uq);
hq=quiver(xq, yq, real(uq), imag(uq), 'k');
figure(2)
hs=surf(xx,yy,rho); shading interp; colormap(jet(256));
axis square; xlim([-1.1 1.1]); ylim([-1.1 1.1]);


dt=0.01/prod(N);
nframes=10000;
for i=1:nframes
    % Update state
    w=solveRK4(w, dt);
    w([1 end],:,:)=w0([1 end],:,:);
    w(:,[1 end],:)=w0(:,[1 end],:);
    u=(w(:,:,2)+1i*w(:,:,3))./w(:,:,1);
    
    % Plot
    uq=interp2(xx',yy',u.',xq,yq); uq=uq./abs(uq);
    hq.UData=real(uq);
    hq.VData=imag(uq);
    hs.ZData=w(:,:,1);
    drawnow;
end

end

function [Fix,Fiy]=iflux(rho, u, v, E)
% Inviscid flux vector
global gamma;
p=(gamma-1)*rho.*(E-(u.^2+v.^2)/2);
H=E+p./rho;

Fix=zeros([size(u), 4]);
Fix(:,:,1)=rho.*u;
Fix(:,:,2)=rho.*u.^2+p;
Fix(:,:,3)=rho.*u.*v;
Fix(:,:,4)=rho.*u.*H;

Fiy=zeros([size(u), 4]);
Fiy(:,:,1)=rho.*v;
Fiy(:,:,2)=rho.*u.*v;
Fiy(:,:,3)=rho.*v.^2+p;
Fiy(:,:,4)=rho.*v.*H;
end

function [Fvx,Fvy]=vflux(u, v)
% Viscous flux vector
global Dx Dy mu vel;
ux=Dx*u;
uy=u*Dy';
vx=Dx*v;
vy=v*Dy';

Fvx=zeros([size(u), 4]);
Fvx(:,:,2)=2/3*mu*(2*ux-vy);
Fvx(:,:,3)=mu*(uy+vx);
Fvx(:,:,4)=Fvx(:,:,2).*u+Fvx(:,:,3).*v+real(vel);

Fvy=zeros([size(u), 4]);
Fvy(:,:,2)=mu*(uy+vx);
Fvy(:,:,3)=2/3*mu*(2*vy-ux);
Fvy(:,:,4)=Fvy(:,:,2).*u+Fvy(:,:,3).*v+imag(vel);
end

function wdot=timeD(w)
global Dx Dy;
rho=w(:,:,1);
u=w(:,:,2)./rho;
v=w(:,:,3)./rho;
E=w(:,:,4)./rho;
[Fix,Fiy]=iflux(rho, u, v, E);
[Fvx,Fvy]=vflux(u, v);
Fx=Fvx-Fix;
Fy=Fvy-Fiy;

wdot=zeros(size(w));
wdot(:,:,1)=Dx*Fx(:,:,1)+Fy(:,:,1)*Dy';
wdot(:,:,2)=Dx*Fx(:,:,2)+Fy(:,:,2)*Dy';
wdot(:,:,3)=Dx*Fx(:,:,3)+Fy(:,:,3)*Dy';
wdot(:,:,4)=Dx*Fx(:,:,4)+Fy(:,:,4)*Dy';
end

function [w]=solveRK4(w, dt)
k1=dt*timeD(w     );
k2=dt*timeD(w+k1/2);
k3=dt*timeD(w+k2/2);
k4=dt*timeD(w+k3  );
w=w+(k1+2*k2+2*k3+k4)/6;
end

function [u]=solenoid(u)
global Dx Dxx Dy Dyy;
divu=Dx*real(u)+imag(u)*Dy';
phi=zeros(size(u));
phi(2:end-1,2:end-1)=sylvester(Dxx(2:end-1,2:end-1), Dyy(2:end-1,2:end-1)', -divu(2:end-1,2:end-1));
u=u+(Dx*phi)+1i*(phi*Dy');
end



