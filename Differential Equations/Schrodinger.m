function u = Schrodinger(N, M)
% Solves the Schrodinger equation in 2D polar coordinates
rho=chebGrid(N);
phi=2*pi*(0:M-1)/M;
[pp, rr]=meshgrid(phi, rho);
[xx, yy]=pol2cart(pp, rr);
u=exp(-10*(abs(rr)-0.4).^2);
dt=6/(N*M);
n=N/2;
figure(1);
h=surf(xx(n:end,:), yy(n:end,:), real(u(n:end,:)), 'EdgeColor', 'none');
shading interp
colormap(jet);
nframes=10000;
for i=1:nframes
    u=solveRK4(u, dt);
    u([1 end],:)=0;
    if(mod(i,10)==1)
        set(h, 'ZData', real(u(n:end,:)));
        drawnow;
    end
end
end

function lap=polarDel2(u)
N=size(u, 1)-1;
Dr=1i*[0:N-1, 0, 1-N:-1]';
th=(1:N-1)'*pi/N; c=cos(th); s=sin(th);
v_hat=fft([u; flip(u(2:N,:), 1)], [], 1);
W1=ifft(bsxfun(@times, Dr, v_hat), [], 1);
W2=ifft(bsxfun(@times, Dr.^2, v_hat), [], 1);
w1=zeros(size(u)); w2=zeros(size(u));
w1(2:N,:)=bsxfun(@times, -1./s, W1(2:N,:));
w2(2:N,:)=bsxfun(@times, s.^-2, W2(2:N,:))-bsxfun(@times, c.*s.^-3, W1(2:N,:));

M=size(u, 2);
Dp=1i*[0:M/2-1, 0, -M/2+1:-1];
u_hat=fft(u, [], 2);
u_phi=ifft(bsxfun(@times, Dp.^2, u_hat), [], 2);

lap=w2+bsxfun(@times, c.^-1, w1)+bsxfun(@times, c.^-2, u_phi);
end

function v=partialTime(u)
v=0.01*polarDel2(u);
end

function u=solveRK4(u, dt)
% Time-stepping by Runge Kutta 4th order.
k1=dt*partialTime(u);
k2=dt*partialTime(u+k1/2);
k3=dt*partialTime(u+k2/2);
k4=dt*partialTime(u+k3);
u=u+(k1+2*k2+2*k3+k4)/6;
end