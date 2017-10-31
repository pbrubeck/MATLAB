function [] = deformBL( n, nframe )
% Example of fluid element deformation on a Boundary Layer

% Select velocity field
delta=@(x) sqrt(x);
P=@(x) 2*x-2*x.^3+x.^4;
vel=@(z) (imag(z)<=delta(real(z))).*P(imag(z)./delta(real(z)))+(imag(z)>delta(real(z)));

% Set grid
x=linspace(0,2,n);
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
xlim([0,4.01]);
ylim([0,4.01]);

i=round(1/5*n)-1;
j=round(1/10*n)-1;
id1=[n*i+j,n*(i+1)+j,n*(i+1)+j+1,i*n+j+1]+1;
h2=patch(real(zz(id1)),imag(zz(id1)),'red');

i=round(4/5*n)-1;
j=round(1/10*n)-1;
id2=[n*i+j,n*(i+1)+j,n*(i+1)+j+1,i*n+j+1]+1;
h3=patch(real(zz(id1)),imag(zz(id1)),'blue');
hold off;

set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
set(gca,'XTick',0:5);
set(gca,'YTick',0:5);
xlabel('$x/L$');
ylabel('$y/\delta(L)$');

view(2);
dt=0.001;
T=2;
tframe=T/nframe;
nframe=0; cp=tframe;
for t=0:dt:T
    zz=solveRK4(zz,dt);
    set(h1,'XData',real(zz));
    set(h1,'YData',imag(zz));
    
    set(h2,'XData',real(zz(id1)));
    set(h2,'YData',imag(zz(id1)));
    set(h3,'XData',real(zz(id2)));
    set(h3,'YData',imag(zz(id2)));
    drawnow;
    if(abs(t-cp)<=dt/2)
        cp=cp+tframe;
        nframe=nframe+1;
        print('-depsc',sprintf('deform%02d_%02d',1,nframe))
    end
end
end