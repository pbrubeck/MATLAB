function [] = schrodCayley(N)
% Solves the Schrodinger equation via Cayley transformations
R1=1;
R2=2;
[D,x]=chebD(N);
r=(R2+R1)/2+(R2-R1)/2*x;
D=2/(R2-R1)*D;
Drr=(diag(r)*D)^2;
[Dtt,th]=periodicD2(N);

Psi=(-(r-R1).*(r-R2).*(r-(R1+R2)/2))*sin(3*th.^2/(2*pi));

[x,w]=ClenshawCurtis(R1,R2,N-1);
W=diag(w);
Psi=Psi/sqrt(sum(r'*W*abs(Psi).^2)*(2*pi/N));


xx=r*cos([0,th]);
yy=r*sin([0,th]);
zz=zeros(N,N+1);

Drr=Drr(2:end-1,2:end-1);
r=r(2:end-1);
V=200*((r-R1)/(R2-R1)).^2;
Hr=-Drr+diag(V.*r.^2);
Ht=-Dtt+diag(30*(th).^2);

frames=2000;
tmax=1;
dt=tmax/frames;
A=diag(r.^2)+1i*dt*Hr;
B=eye(size(Ht))+1i*dt*Ht';
AA=conj(A);
BB=conj(B);

W=W(2:end-1,2:end-1);
h=surf(xx,yy,zz);
shading interp; colormap(jet(256));
zlim([0, 0.6]); view(0,90);

[Q,U]=eig(A,'vector');
[P,V]=eig(B,'vector');
for t=0:dt:tmax
    C=AA*Psi(2:end-1,:)+Psi(2:end-1,:)*BB;
    Psi(2:end-1,:)=solve(U,Q,V,P,C);
    zz=abs(Psi).^2;
    set(h,'ZData',[zz(:,end), zz]);
    drawnow;
    %disp(sum(r'*W*abs(Psi(2:end-1,:)).^2)*(2*pi/N));
end

end

function [X]=solve(U,Q,V,P,C)
Z=zeros(size(C));
G=Q\C*P;
n=size(P,1);
for i=n:-1:1
    Z(:,i)=G(:,i)./(U+V(i));
end
X=Q*Z/P;
end

