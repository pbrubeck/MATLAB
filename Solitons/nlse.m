function [] = nlse(N, t1, t2, L, init, method)
% NonLinear Schrodinger Equation in one dimension
% 1i u_t + 1/2 u_xx +|u|^2 u = 0
%
% N: number of collocation points
% t1: initial time
% t2: final time
% init: bright, dark1, dark2, peregrine, breather
% method: cheb, fft, herm

nframes=1024;
dt=1/1024;
m=ceil((t2-t1)/(dt*(nframes-1)));
dt=(t2-t1)/(m*(nframes-1));
Amax=2;

% Linear propagator, half step
if nargin<6, method='cheb'; end;
switch(method)
    case 'cheb', [x,w,D,Q]=nlsecheb(N,dt/2);
    case 'herm', [x,w,D,Q]=nlseherm(N,dt/2);
    case 'fft',  [x,w,D,Q]=nlsefft(N, dt/2);
end

beta=-1;
switch(init)
    case {0, 'bright'}
        u=nlsebright(2,-5,0,x,t1)+nlsebright(-2,5,pi/3,x,t1);
    case {1, 'dark1'}
        beta=1;
        u=nlsedark(0, pi/4, -5, 0, x, t1)+nlsedark(0, -pi/4, 5, 0, x, t1);
    case {2, 'dark2'}
        beta=1;
        k=0.2;
        u=nlsedark(k, pi/4, -5, 0, x, t1)+nlsedark(k, -pi/4, 5, 0, x, t1);
    case {3, 'peregrine'}
        u=peregrine(x, t1);
    case {4, 'breather'}
        u=nlsebreather(4, x, t1);
end

Z=w'*(abs(u).^2);                                        display(Z);
P=1i/2*w'*(u.*linprop(D,conj(u))-conj(u).*linprop(D,u)); display(P);
E=w'*(abs(linprop(D,u)).^2+abs(u).^4);                   display(E);

f1=figure(1);
h=plot(x, abs(u), 'b', 'LineWidth', 2);
xlim([-L,L]); ylim([0,Amax]); axis manual;
set(gca, 'XTick', [-L,0,L]);
set(gca, 'YTick', [0,Amax/2,Amax]);
xlabel('x'); title('|\Psi|');
set(gca,'FontSize',36);
drawnow;

% Time propagation
inc=0; tinc=0; tframe=(t2-t1)/8;
udata=zeros(nframes,length(u));
udata(1,:)=u;
for i=2:nframes
    for j=1:m
        u=linprop(Q,u);
        u=u.*exp(-1i*beta*dt*(abs(u).^2));
        u=linprop(Q,u);
        
        tinc=tinc+dt;
        if(tinc>=tframe)
            print(f1,sprintf('%s\\data\\%s%02d',pwd,init,inc),'-depsc');
            tinc=0; inc=inc+1;
        end
    end
    udata(i,:)=u;
    set(h, 'YData', abs(u));
    drawnow;
end
print(f1,sprintf('%s\\data\\%s%02d',pwd,init,inc),'-depsc');

Z=w'*(abs(u).^2);                                        display(Z);
P=1i/2*w'*(u.*linprop(D,conj(u))-conj(u).*linprop(D,u)); display(P);
E=w'*(abs(linprop(D,u)).^2+abs(u).^4);                   display(E);

t=t1:m*dt:t2;

figure(2);
roi=abs(x)<=L;
surf(x(roi),t,abs(udata(:,roi)).^2);
colormap(jet(256)); colorbar(); shading interp; view(2);
xlabel('x'); ylabel('t'); title('|\Psi|^2');

figure(3);
imagesc(x(roi),t,angle(udata(:,roi)));
colormap(hsv(256)); colorbar(); shading interp;
xlabel('x'); ylabel('t'); title('arg(\Psi)');
end

function u=linprop(Q,u)
if isfloat(Q)
    u=Q*u;
else
    u=Q(u);
end
end