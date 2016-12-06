function [] = heatExp( N )
dt=0.001;

[x,w]=gauleg(-1,1,N); x=x(:); w=w(:);
V=exp(1i*pi/2*x*(-N/2:N/2-1))/sqrt(2);
L=-pi^2/4*(-N/2:N/2-1)'.^2;

u=x-x.^2;
uhat=V'*(w.*u);

figure(1);
h=plot(x, real(u)); axis manual;

nframes=1000;
for i=1:nframes
    u=V*(exp(i*dt*L).*uhat);
    set(h, 'YData', real(u));
    drawnow;
end
end