function u=NavierStokes(n)
% Solves for the velocity field assuming periodic boundary conditions.
nu=1.00;
dt=0.0001;

% Construct 3D grid
gv=2*pi*(0:n-1)/n;
[xx,yy,zz]=meshgrid(gv);

% Initial velocity field
u=cat(4, -sin(xx-pi), -cos(yy-pi), -sin(zz-pi));

% Initialize spectral differential operators
ii=[0:n/2-1, -n/2:-1]';
jj=[0:n/2-1, -n/2:-1]';
kk=[0:n/2-1, -n/2:-1]';
Dx=1i*ii;
Dy=1i*reshape(jj, [1 n]);
Dz=1i*reshape(kk, [1 1 n]);

figure(1); % Initialize quiver plot
h=quiver3(xx, yy, zz, u(:,:,:,1), u(:,:,:,2), u(:,:,:,3));
axis equal; view(3);

nframes=10000;
for t=1:nframes
    tic
    u=solveRK4(u, Dx, Dy, Dz, nu, dt); 
    title(sprintf('%.0f fps', 1./toc));
    if(mod(t,10)==0)
        set(h, 'UData', u(:,:,:,1));
        set(h, 'VData', u(:,:,:,2));
        set(h, 'WData', u(:,:,:,3));
        drawnow;
    end
end
end

function ut=partialTime(u, Dx, Dy, Dz, nu)
u1=u(:,:,:,1);    u2=u(:,:,:,2);    u3=u(:,:,:,3);
v1=fft(u, [], 1); v2=fft(u, [], 2); v3=fft(u, [], 3);
% Jacobian tensor = parital derivatives on each direction
Jx=ifft(bsxfun(@times, Dx, v1), [], 1, 'symmetric');
Jy=ifft(bsxfun(@times, Dy, v2), [], 2, 'symmetric');
Jz=ifft(bsxfun(@times, Dz, v3), [], 3, 'symmetric');
% Advection = Matrix product of Jacobian with the field
adv=bsxfun(@times, u1, Jx)+bsxfun(@times, u2, Jy)+bsxfun(@times, u3, Jz);
% Laplacian = Sum of second partial derivatives
uxx=ifft(bsxfun(@times, nu*Dx.^2, v1), [], 1, 'symmetric');
uyy=ifft(bsxfun(@times, nu*Dy.^2, v2), [], 2, 'symmetric');
uzz=ifft(bsxfun(@times, nu*Dz.^2, v3), [], 3, 'symmetric');
lap=uxx+uyy+uzz;
ut=lap-adv;
end

function u=solveRK4(u, Dx, Dy, Dz, nu, dt)
% Time-stepping by Runge Kutta 4th order
k1=dt*partialTime(u,      Dx, Dy, Dz, nu);
k2=dt*partialTime(u+k1/2, Dx, Dy, Dz, nu);
k3=dt*partialTime(u+k2/2, Dx, Dy, Dz, nu);
k4=dt*partialTime(u+k3,   Dx, Dy, Dz, nu);
u=u+(k1+2*k2+2*k3+k4)/6;
end