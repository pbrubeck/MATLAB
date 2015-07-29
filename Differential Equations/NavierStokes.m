function u=NavierStokes(n)
% Solves for the velocity field assuming periodic boundary conditions.

nu=0.75;
dt=0.001;

i=[0:n/2-1, 0, -n/2+1:-1];
j=[0:n/2-1, 0, -n/2+1:-1];
k=[0:n/2-1, 0, -n/2+1:-1];
[Dx,Dy,Dz]=meshgrid(1i*i, 1i*j, 1i*k);
[i2,j2,k2]=meshgrid(i.^2, j.^2, k.^2);
D2=-nu*(i2+j2+k2);

x=2*pi*(0:n-1)/n;
y=2*pi*(0:n-1)/n;
z=2*pi*(0:n-1)/n;
[xx,yy,zz]=meshgrid(x, y, z);

u=cat(4, sin(xx), cos(yy).*sin(zz), sin(yy));

figure(1);
h=quiver3(xx, yy, zz, u(:,:,:,1), u(:,:,:,2), u(:,:,:,3));
axis equal; view(-26, 32);

nframes=1000;
for t=1:nframes
    tic
    u=solveRK4(u, dt, Dx, Dy, Dz, D2);
    title(sprintf('Calculation time %.0f ms', 1000*toc));
    if(mod(t,10)==0)
        set(h, 'UData', u(:,:,:,1));
        set(h, 'VData', u(:,:,:,2));
        set(h, 'WData', u(:,:,:,3));
        drawnow;
    end
end
end

function ut=partialTime(u, Dx, Dy, Dz, D2)
u_hat=fftn(u(:,:,:,1));
v_hat=fftn(u(:,:,:,2));
w_hat=fftn(u(:,:,:,3));

% Laplacian
lap=cat(4, ifftn(D2.*u_hat), ifftn(D2.*v_hat), ifftn(D2.*w_hat));

% Jacobian
Jx=cat(4, ifftn(Dx.*u_hat), ifftn(Dx.*v_hat), ifftn(Dx.*w_hat));
Jy=cat(4, ifftn(Dy.*u_hat), ifftn(Dy.*v_hat), ifftn(Dy.*w_hat));
Jz=cat(4, ifftn(Dz.*u_hat), ifftn(Dz.*v_hat), ifftn(Dz.*w_hat));

% Advection
adv  =  bsxfun(@times, u(:,:,:,1), Jx);
adv=adv+bsxfun(@times, u(:,:,:,2), Jy);
adv=adv+bsxfun(@times, u(:,:,:,3), Jz);
ut=lap-adv;
end

function u=solveRK4(u, dt, Dx, Dy, Dz, D2)
k1=dt*partialTime(u,      Dx, Dy, Dz, D2);
k2=dt*partialTime(u+k1/2, Dx, Dy, Dz, D2);
k3=dt*partialTime(u+k2/2, Dx, Dy, Dz, D2);
k4=dt*partialTime(u+k3,   Dx, Dy, Dz, D2);
u=u+(k1+2*k2+2*k3+k4)/6;
end