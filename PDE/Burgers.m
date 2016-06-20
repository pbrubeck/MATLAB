function [] = Burgers(N)
% p27.m - Solve KdV eq. u_t + uu_x + u_xxx = 0 on [-pi,pi] by
%         FFT with integrating factor v = exp(-ik^3t)*u-hat.

% Set up grid and two-soliton initial data:
nu = 0.01;
dt = .4/N^2;
x = (2*pi/N)*(-N/2:N/2-1)';
u = zeros(size(x));
u(x<0)=sin(x(x<0)).^2;
v = fft(u);
k = [0:N/2-1 0 -N/2+1:-1]'; k2 = k.^2;
E = exp(-nu*dt*k2/2); E2 = E.^2;

% Solve PDE and plot results:
tmax = 10; 
nplt = floor((tmax/250)/dt); 
nmax = round(tmax/dt);

h=plot(x, u);
xlim([-pi,pi]);
ylim([0,1]);
for n = 1:nmax
    t = n*dt; 
    g = -.5i*dt*k;
    a = g.*fft(real( ifft(     v    ) ).^2);
    b = g.*fft(real( ifft(E.*(v+a/2)) ).^2)./E;     
    c = g.*fft(real( ifft(E.*(v+b/2)) ).^2)./E;
    d = g.*fft(real( ifft(E2.*(v+c)) ).^2)./E2;
    v = E2.*(v+(a+2*b+2*c+d)/6);
    if mod(n,nplt) == 0 
          u = real(ifft(v));
          set(h, 'YData', u);
          drawnow;
    end
end

end