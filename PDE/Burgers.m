function [] = Burgers(N)
% Solve Burgers eq. u_t + uu_x - nu*u_xx = 0 on [-pi,pi] by
% FFT with integrating factor v = exp(-nu*k^2*t)*u-hat.

% Set up grid and two-soliton initial data:
nu=0.01;
dt=0.4/N^2;
x=(2*pi/N)*(-N/2:N/2-1)';
u=zeros(size(x));
u(x<0)=sin(x(x<0)).^2;
v=fft(u);
k=[0:N/2-1 0 -N/2+1:-1]';
E=exp(-nu*dt*k.^2/2); E2 = E.^2;
g=-0.5i*dt*k;

% Solve PDE and plot results:
tmax=20; 
nplt=floor((tmax/1000)/dt); 
nmax=round(tmax/dt);

figure(1);
h=plot(x, u);
xlim([-pi,pi]);
ylim([0,1]);
for n = 1:nmax 
    k1=g.*fft(real(ifft(v)).^2);
    k2=g.*fft(real(ifft(E.*(v+k1/2))).^2)./E;     
    k3=g.*fft(real(ifft(E.*(v+k2/2))).^2)./E;
    k4=g.*fft(real(ifft(E2.*(v+k3))).^2)./E2;
    v = E2.*(v+(k1+2*k2+2*k3+k4)/6);
    if mod(n,nplt) == 0 
          u = real(ifft(v));
          set(h, 'YData', u);
          drawnow;
    end
end

end