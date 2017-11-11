function [] = wavePML1( N )
% Solves the wave equation in 1D using a Perfectly Matched Layer for BCs.
% Uses the method of lines for time evolution with Chebyshev
% discretization.

% Differential operator
[D,x]=chebD(N);
dt=1e-2;

% Layer
xl=0.95;
layer=abs(x)>xl;
sigma=zeros(N,1);
sigma(layer)=N^(3/2)*((abs(x(layer))-xl)/(1-xl)).^3;

% Linear propagator
A1=-diag(sigma(2:end-1));
A2=D(2:end-1,2:end-1);
A=[A1, A2; A2, A1];
Q=expm(A*dt);

% Initial conditions
a=30;
x0=-0.1;
u0=1./(1+(a*(x-x0)).^2);
v0=-u0;
w=[u0(2:end-1), v0(2:end-1)];

% Plot
figure(1);
h=plot(x, u0, 'b', 'Linewidth', 1);
xlim([-xl,xl]); ylim([-1,1]); xlabel('x'); title('\Psi(x)');
nframes=2*ceil((xl-x0)/dt)+1;
wdata=zeros(N,nframes);
for i=1:nframes
    wdata(:,i)=[0; w(:,1); 0];
    w=reshape(Q*w(:), size(w));
    set(h, 'YData', [0; w(:,1); 0]);
    drawnow;
end
t=dt*(0:(nframes-1));
[xx,tt]=ndgrid(x,t);
figure(2);
surf(xx,tt,wdata); shading interp;  colormap(jet(256)); view(2); 
xlim([-xl,xl]); ylim([0,dt*(nframes-1)]);
xlabel('x'); ylabel('t'); title('\Psi(x,t)');
end