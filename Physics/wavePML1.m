function [] = wavePML1( N )
% Solves the wave equation in 1D using a Perfectly Matched Layer for BCs.
% Uses the method of lines for time evolution with Chebyshev
% discretization.
[D,x]=chebD(N);
dt=5e-2;

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
x0=0.5;
u0=exp(-(10*(x-x0)).^2)-exp(-(10*(x+x0)).^2);
v0=exp(-(10*(x-x0)).^2)+exp(-(10*(x+x0)).^2);
w=[u0(2:end-1), v0(2:end-1)];

% Plot
figure(1);
h=plot(x, u0, 'b', 'Linewidth', 1);
xlim([-xl,xl]); ylim([-1,1]);
nframes=ceil((xl+x0)/dt);
for i=1:nframes
    w=reshape(Q*w(:), size(w));
    set(h, 'YData', [0; w(:,1); 0]);
    drawnow;
    if(max(w(1,:))>1)
        pause();
    end
end
end