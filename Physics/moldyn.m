function [Ti,Ei,ti] = moldyn(N,d)
% N number of particles
% d number of dimensions

filename = sprintf('moldyn%d.gif', N);

L=4*N^(1/d);    % Lenght of container
T0=0;           % Velocity variance
R=2^(1/6)/2;    % Collision radius (Fixed)
dL=R/1e5;       % Max step
tf=10;          % Final time
fr=10;          % Framerate
nframes=1+ceil(tf*fr);

% Initial configuration
E=1;
while E>=0
    x=L*rand(N,d);
    v=sqrt(T0)*randn(N,d);
    [U,F]=lennardjones(x,box(x,L),L);
    E=sum(U)/2+sum(sum(v.^2))/2;
    T=sum(sum(v.^2))/(N*d);
end

% Plot spheres
[s1,s2,s3]=sphere(20);
[s1,s2,s3]=deal(R*s1,R*s2,R*s3);
cmap=hsv(N);
figure(1); 
h=cell(N,1);
for j=1:N
    h{j}=surf(s1+x(j,1), s2+x(j,2), s3+x(j,3), ...
        'EdgeColor', 'none', 'FaceColor', cmap(j,:));
    if(j==1), hold on; end;
end
camlight; hold off;
xlim([0,L]); ylim([0,L]); zlim([0,L]);
daspect([1,1,1]);

t=0;
title(sprintf('T=%.2f, E=%.2f, t=%.2f',T,E,t));
print('hw08g01','-dpng');

% Create .gif file
im=frame2im(getframe(gcf));
[imind,cm]=rgb2ind(im,256);
imwrite(imind,cm,filename,'gif','DelayTime',0,'Loopcount',inf);

ti=zeros(nframes,1); ti(1)=t;
Ei=zeros(nframes,1); Ei(1)=E;
Ti=zeros(nframes,1); Ti(1)=T;

for i=2:nframes
    while(t<(i-1)/fr)
        % Compute time step
        dt=min(1/fr, sqrt(dL/max(sqrt(dot(F,F,2)))));
        % Velocity Verlet method
        v=v+F*dt/2;
        x=mod(x+v*dt,L);
        [U,F]=lennardjones(x,box(x,L),L);
        v=v+F*dt/2;
        t=t+dt;
    end
    
    K=sum(sum(v.^2))/2;
    E=K+sum(U)/2;
    T=2*K/(N*d);
    Ti(i)=T; Ei(i)=E; ti(i)=t;
    title(sprintf('T=%.2f, E=%.2f, t=%.2f',T,E,t));
    for j=1:N
        set(h{j}, 'XData', s1+x(j,1));
        set(h{j}, 'YData', s2+x(j,2));
        set(h{j}, 'ZData', s3+x(j,3));
    end
    drawnow;

    im=frame2im(getframe(gcf));
    [imind,cm]=rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','DelayTime',0,'WriteMode','append');
end

print('hw08g02','-dpng');
plot(ti,Ti,'Linewidth',2);
title('Average temperature');
xlabel('time');
print('hw08g03','-depsc');
end

function x1=box(x0,L)
[dx,dy,dz]=ndgrid([-L,0,L]);
A=kron([dx(:),dy(:),dz(:)],ones(size(x0,1),1));
x1=A+repmat(x0,[3^size(x0,2),1]);
end

function [U,F]=lennardjones(x0,x1,L)
EPS=10;
F=zeros(size(x0));
dx=zeros(size(x0,1),size(x1,1),size(x0,2));
for i=1:size(x0,2)
    [xx0,xx1]=ndgrid(x0(:,i),x1(:,i));
    dx(:,:,i)=xx0-xx1;
end
r2=sum(dx.^2, 3);
r6=r2.^3; r12=r6.^2;
U=4*EPS*(1./r12-1./r6);
A=4*EPS*(12./r12-6./r6)./r2;
B=(r2>L^2/4)|(r2==0);
A(B)=0;
U(B)=0; U=sum(U,2);
for i=1:size(x0,2)
    F(:,i)=dot(A, dx(:,:,i), 2);
end
end