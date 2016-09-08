function [] = schrodExp( N )
[D,x]=chebD(N);
[~,w]=ClenshawCurtis(-1,1,N-1);

dt=0.001;
V0=1000;
T=-1/2*D*D;
V=diag(V0*x.^2);
H=T+V;

[V,L]=eig(H(2:end-1,2:end-1),'vector');
L=[0;0;L];
U=zeros(N,N);
U(:,1:2)=null(H);
U(2:end-1,3:end)=V;
Q=U*diag(exp(-1i*dt*L))/U;

sigma=0.8;
k=1/sigma;
Psi0=normc(exp(-(x/sigma).^2-1i*k*x), w);
h=plot(x,abs(Psi0).^2); drawnow; ylim([0, 4]);
for i=1:10000
    Psi0=Q*Psi0;
    set(h, 'YData', abs(Psi0).^2);
    drawnow;
end
end

