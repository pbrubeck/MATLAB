function [] = deform( n, mode, nframe )
% Example of fluid element deformation

% Select velocity field
U=@(r) r.^2;
V=@(r) 1./(r.^2+1);
switch(mode)
    case 1
        vel=@(z) U(real(z))+1i*V(imag(z));
    case 2
        vel=@(z) U(imag(z));
    case 3
        vel=@(z) imag(z);
    case 4
        vel=@(z) U(abs(z)).*(1i*z)./abs(z);
    case 5
        vel=@(z) (1i*z);
    case 6
        vel=@(z) (1i*z)./(abs(z).^2);
end

% Set grid
x=linspace(-1,1,n);
[xx,yy]=ndgrid(x);
zz=xx+1i*yy;

% Runge-Kutta
function z=solveRK4(z, dt)
    k1=dt*vel(z);
    k2=dt*vel(z+k1/2);
    k3=dt*vel(z+k2/2);
    k4=dt*vel(z+k3);
    z=z+(k1+2*k2+2*k3+k4)/6;
end

figure(1); clf; hold on;
h1=mesh(real(zz),imag(zz),0*real(zz)-1,'linewidth',1);
colormap([0,0,0]);
xlim([-2.5,2.5]);
ylim([-2.5,2.5]);

i=round(2/3*n)-1;
id=[n*i+i,n*(i+1)+i,n*(i+1)+i+1,i*n+i+1]+1;
h2=patch(real(zz(id)),imag(zz(id)),'red');
hold off;

set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gca,'XTick',-2:2);
set(gca,'YTick',-2:2);
xlabel('$x$');
ylabel('$y$');

view(2);
dt=0.001;
T=2;
tframe=T/nframe;
nframe=0; cp=tframe;
for t=0:dt:T
    zz=solveRK4(zz,dt);
    set(h1,'XData',real(zz));
    set(h1,'YData',imag(zz));
    
    set(h2,'XData',real(zz(id)));
    set(h2,'YData',imag(zz(id)));
    drawnow;
    if(abs(t-cp)<=dt/2)
        cp=cp+tframe;
        nframe=nframe+1;
        %print('-depsc',sprintf('deform%02d_%02d',mode,nframe))
    end
end
end