function [] = waveExp(N)
a=[1,1]; b=[0,0];
[D,x]=chebD(N); D2=D*D;
[xx,yy]=ndgrid(x);
[A1, G1, H1, C1] = setBC(D2, D, a(1), b(1)*[1,-1]);
[A2, G2, H2, C2] = setBC(D2, D, a(2), b(2)*[1,-1]);
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

U0=(xx.^2+yy.^2)/5;
V0=-(0*D*U0+0*U0*D');

kd=2:N-1;
rd=[1,N];

% Solve Laplace
F=-D2(:,rd)*U0(rd,:)-U0(:,rd)*D2(:,rd)';
UH=zeros(N);
UH(kd,kd)=W1(kd,kd)*((W1(kd,kd)\F(kd,kd)/W2(kd,kd)')./LL(kd,kd))*W2(kd,kd)';
UH(rd,:)=H1*U0(rd,:)+G1*UH(kd,:);
UH(:,rd)=U0(:,rd)*H2'+UH(:,kd)*G2';

UU=W1\(U0-UH)/W2';
VV=W1\(V0)/W2';
U=@(t) UH+W1*(UU.*cos(ww*t)+VV.*sinc(ww*t)*t)*W2';

figure(1);
h=surf(xx,yy,U0);
colormap(jet(256));
zlim([-0.1,0.6]); 
axis square; alpha(0.9);
camlight; shading interp; 

dt=2*pi/min(min(ww(ww>0)))/(100);
nframes=101;
for i=0:nframes-1
    set(h,'ZData',U(i*dt));
    drawnow;
end
end