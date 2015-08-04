function u = wavePDE(N)
% Solves the wave equation in 2D
[xx,yy]=meshgrid(chebGrid(N));
u=exp(-40*((xx-.4).^2+yy.^2));
u_old=u;
dt=6/N^2;
h=surf(xx, yy, u, 'EdgeColor', 'none');
colormap(jet);
nframes=10000;
for i=1:nframes
    lap=chebfftD2(u,1)+chebfftD2(u,2);
    u_new=2*u-u_old+dt^2*lap;
    u_old=u; u=u_new;
    u([1 end],:)=0; u(:,[1 end])=0;
    if(mod(i,10)==1)
        set(h, 'ZData', u);
        drawnow;
    end
end
end