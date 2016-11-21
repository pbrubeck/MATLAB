function [] = nlse2(N, init, method)
% NonLinear Schrodinger Equation in two dimensions
% 1i u_t + 1/2 u_xx +|u|^2 u = 0

t0=-pi; tf=pi; 
nframes=1024;
dt=1/1024;
m=ceil((tf-t0)/(dt*(nframes-1)));
dt=(tf-t0)/(m*(nframes-1));
x0=3;

% Linear propagator, half step
if nargin==2, method='cheb'; end;
switch(method)
    case 'cheb', [Q,x]=nlsecheb(N,dt/2);
    case 'herm', [Q,x]=nlseherm(N,dt/2);
    case 'fft',  [Q,x]=nlsefft(N, dt/2);
end
[xx,yy]=ndgrid(x);

if(init==0)
    U=peregrine(hypot(xx,yy),t0);
else
    U=kmbreather(init,hypot(xx,yy),t0);
end

figure(1);
h=surf(xx, yy, abs(U).^2);
shading interp; colormap(jet(256));
xlim([-x0,x0]); ylim([-x0,x0]); zlim([0,10]); axis manual;
xlabel('x'); ylabel('y'); title('|\Psi|^2');
drawnow;

for i=2:nframes
    for j=1:m
        U=linprop2(Q,U);
        U=U.*exp(1i*dt*(abs(U).^2));
        U=linprop2(Q,U);
    end
    set(h, 'ZData', abs(U).^2);
    drawnow;
end

end

function U=linprop2(Q,U)
if isfloat(Q)
    U=Q*U*Q.';
else
    U=Q(Q(U).').';
end
end