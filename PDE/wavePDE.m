function u = wavePDE(N)
% Solves the wave equation in 2D.
% Example on how to use Chebyshev spectral methods in two spatial dimensions
% and classic Runge-Kutta for time evolution of second order PDE.
[xx,yy]=ndgrid(chebGrid(N));
u(:,:,1)=exp(-40*((xx-.4).^2+yy.^2));
u(:,:,2)=0;
dt=6/N^2;
h=surf(xx, yy, u(:,:,1), 'EdgeColor', 'none');
colormap(jet(256)); alpha(0.85); shading interp;
view(2); zlim([-1,1]); axis square; 

nframes=10000;
for i=1:nframes
    u=solveRK4(u,dt);
    u([1 end],:,:)=0; u(:,[1 end],:)=0;
    if(mod(i,2)==1)
        set(h, 'ZData', real(u(:,:,1)));
        drawnow;
    end
end
end

function v=partialTime(u)
% Returns the vector with the first and second temporal derivatives.
v(:,:,1)=u(:,:,2);
v(:,:,2)=chebfftD2(u(:,:,1),1)+chebfftD2(u(:,:,1),2);
end

function u=solveRK4(u, dt)
% Time-stepping by Runge Kutta 4th order.
k1=dt*partialTime(u);
k2=dt*partialTime(u+k1/2);
k3=dt*partialTime(u+k2/2);
k4=dt*partialTime(u+k3);
u=u+(k1+2*k2+2*k3+k4)/6;
end