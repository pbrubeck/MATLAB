function [] = waveExp(N)
a=[1,1]; b=[0,0];
[D,x]=chebD(N); D2=D*D;
[xx,yy]=ndgrid(x);
[A1, G1] = setBC(D2, D, a(1), b(1)*[1,-1]);
[A2, G2] = setBC(D2, D, a(2), b(2)*[1,-1]);
W1=zeros(N); L1=zeros(N,1);
W2=zeros(N); L2=zeros(N,1);
W1(:,[1,end])=null(D2);
W2(:,[1,end])=null(D2);
[W1(2:end-1,2:end-1), L1(2:end-1)] = eig(A1, 'vector');
[W2(2:end-1,2:end-1), L2(2:end-1)] = eig(A2, 'vector');
W1([1,end],2:end-1)=G1*W1(2:end-1,2:end-1);
W2([1,end],2:end-1)=G2*W2(2:end-1,2:end-1);
[L1, L2] = ndgrid(L1, L2);
LL=L1+L2; ww=sqrt(-LL);
ww(L1==0)=0;
ww(L2==0)=0;

zz=xx+1i*yy;
U0=real(zz.^2)+1*exp(-50*(xx.^2+yy.^2)/2);
V0=-(0*D*U0-1*U0*D');

UU=W1\U0/W2';
VV=(W1\V0/W2')./ww;
VV(ww==0)=0;

U=@(t) W1*(UU.*cos(ww*t)+VV.*sin(ww*t))*W2';

figure(1);
h=surf(xx,yy,U0);
colormap(jet(256));
axis square manual;
camlight; shading interp; 

dt=2*pi/min(min(ww(ww>0)))/(100);
nframes=1001;
for i=0:nframes-1
    set(h,'ZData',U(i*dt));
    drawnow;
end
end