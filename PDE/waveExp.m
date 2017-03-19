function [] = waveExp(N)
a=[1,1;1,1]; 
b=0*[1,1;1,1]; 

rd=[1,N];
kd=2:N-1;
[x,w]=ClenshawCurtis(-1,1,N); w=w(:);
[D,x]=chebD(N); D2=D*D;
[xx,yy]=ndgrid(x);
[A1,G1]=setBC(D2, D, a(1,:), b(1,:));
[A2,G2]=setBC(D2, D, a(2,:), b(2,:));
W1=zeros(N); L1=zeros(N,1);
W2=zeros(N); L2=zeros(N,1);
W1(:,rd)=null(D2);
W2(:,rd)=null(D2);
[W1(kd,kd), L1(kd)] = eig(A1, 'vector');
[W2(kd,kd), L2(kd)] = eig(A2, 'vector');
W1(rd,kd)=G1*W1(kd,kd);
W2(rd,kd)=G2*W2(kd,kd);
[L1, L2] = ndgrid(L1, L2);
LL=L1+L2; ww=sqrt(-LL);
ww(rd,:)=0;
ww(:,rd)=0;

zz=xx+1i*yy;
U0=real(zz.^2)+1*exp(-100*(xx.^2+yy.^2)/2);
V0=-(1*D*U0+0*U0*D');

UU=pinv(W1)*U0*pinv(W2');
VV=(pinv(W1)*V0*pinv(W2'))./ww;
VV(rd,:)=0;
VV(:,rd)=0;

U=@(t) W1*(UU.*cos(ww*t)+VV.*sin(ww*t))*W2';
V=@(t) W1*(ww.*(UU.*sin(ww*t)-VV.*cos(ww*t)))*W2';

figure(1);
h=surf(xx,yy,U0);
colormap(jet(256));
zlim([-2,2]);
axis square manual;
camlight; shading interp; 

dt=2*pi/min(min(ww(ww>0)))/(100);
nframes=1001;
for i=0:nframes-1
    uu=U(i*dt);
    vv=V(i*dt);
    set(h,'ZData',uu);
    E=1/2*w'*(vv.^2+(D*uu).^2+(uu*D').^2)*w;
    title(sprintf('E=%.10f',E));
    drawnow;
end
end