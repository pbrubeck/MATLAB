function [] = pwePML1( N )
% Solves the paraxial wave equation in 1D using a Perfectly Matched Layer for BCs.
% Uses the method of lines for time evolution with Hermite discretization.
lambda=100;
k=2*pi/lambda;
[D,x]=hermD(N);
dt=5e-4;

% Layer
xl=3/4*x(end);
layer=abs(x)>xl;
sig=zeros(N,1);
sig(layer)=N^(3/2)*((abs(x(layer))-xl)/(x(end)-xl)).^3;

% Linear propagator
D2=D*D;
A11=1i/(2*k)*D2(2:end-1,2:end-1);
A12=-(1i/(2*k)*D(2:end-1,2:end-1)+eye(N-2))*diag(sig(2:end-1));
A21=D(2:end-1,2:end-1);
A22=-diag(sig(2:end-1));
A=[A11, A12; A21, A22];
Q=expm(A*dt);

% Initial conditions
u0=cos(2*pi*x).*exp(-x.^2/2);
v0=zeros(size(x));
w=[u0(2:end-1),v0(2:end-1)];

% Plot
figure(1);
h=plot(x, u0, 'b', 'Linewidth', 1);
xlim([-xl,xl]); ylim([-1,1]);
nframes=ceil(xl/(lambda*dt));
for i=1:nframes
    w=reshape(Q*w(:), size(w));
    set(h, 'YData', [0; abs(w(:,1)).^2; 0]);
    drawnow;
end
end