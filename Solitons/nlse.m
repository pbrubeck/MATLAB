function [] = nlse(init,Q,x,dt,t0,tf,p)
% NLSE NonLinear Schrodinger Equation
% 1i u_t + 1/2 u_xx +|u|^2 u = 0
% T=-1/2*(D*D); Q=expm(-1i*dt/2*T);
if(nargin<7)
    p=1;
end

if(init==0)
    psi=peregrine(x,t0);
else
    psi=kmbreather(init,x,t0);
end

nframes=1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
x0=3;

figure(1);
h=plot(x, abs(psi), 'b', 'LineWidth', 2);
xlim([-x0,x0]); ylim([0,4]); axis manual;
xlabel('x'); title('|\Psi|');
drawnow;

udata=zeros(nframes,length(psi));
udata(1,:)=psi;
u=psi./p;
for i=2:nframes
    for j=1:m
        u=Amtimes(Q,u);
        u=u.*exp(1i*dt*(abs(p.*u).^2));
        u=Amtimes(Q,u);
    end
    psi=p.*u;
    udata(i,:)=psi;
    set(h, 'YData', abs(psi));
    drawnow;
end

figure(2);
roi=abs(x)<=x0;
surf(x(roi),t0:m*dt:tf,abs(udata(:,roi)).^2);
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('|\Psi|^2');
end

function x=Amtimes(A,x)
if isfloat(A)
    x=A*x;
else
    x=A(x);
end
end