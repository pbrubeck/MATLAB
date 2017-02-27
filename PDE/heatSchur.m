function [] = heatSchur( N )
[D,x]=chebD(N);
D2=D*D;
xx=[x+1; x(2:end)-1];

s1=0.02;
s2=1;

A11=s1*D2(2:end-1,2:end-1);
A13=[s1*D2(2:end-1,[1,end]), zeros(N-2,1)];

A22=s2*D2(2:end-1,2:end-1);
A23=[zeros(N-2,1), s2*D2(2:end-1,[1,end])];

A31=[zeros(1,N-2); D(end,2:end-1); zeros(1,N-2)];
A32=[zeros(1,N-2); -D(1 ,2:end-1); zeros(1,N-2)];
A33=[1,0,0; D(end,1), D(end,end)-D(1,1), -D(1,end); 0,0,1];
S=A33-A31*(A11\A13)-A32*(A22\A23);

G1=-A33\A31;
G2=-A33\A32;
G=[G1, G2];
A=[A11+A13*G1, A13*G2; A23*G1, A22+A23*G2];

AA=zeros(2*N-1);
AA(1:N-1,1:N)=s1*D2(1:N-1,:);
AA(N+1:end,N:end)=s2*D2(2:N,:);
AA(N,1:N-1)=s1*D2(N,1:N-1); 
AA(N,N)=s1*D2(N,N)-s2*D2(1,1); % Continuity on operator
AA(N,N+1:end)=-s2*D2(1,2:N);

V=zeros(2*N-1);
L=zeros(2*N-1,1);
rd=[1, N, 2*N-1];
kd=[2:N-1, N+1:2*N-2];
[V(kd,kd), L(kd)]=eig(A,'vector');
V(rd,kd)=G*V(kd,kd);
V(:,rd)=null(AA); % nullspace

dt=0.01;
Q=V*diag(exp(dt*L))/V;


u=(exp(-(xx+1).^2)+exp(-(xx-1).^2));

figure(1);
h=plot(xx, u); drawnow; axis manual;

nframes=10000;
for i=1:nframes
    u=Q*u;
    set(h, 'YData', u);
    drawnow;
end
end

