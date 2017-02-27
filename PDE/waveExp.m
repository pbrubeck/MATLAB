function [] = waveExp(N)
c=0;
[D,x]=chebD(N); D2=D*D;
[A1, G1, H1] = setBC(D2, D, 1, c*[1,-1]);
[A2, G2, H2] = setBC(D2, D, 1, c*[1,-1]);
W1=zeros(N); L1=zeros(N,1);
W2=zeros(N); L2=zeros(N,1);
W1(:,[1,end])=null(D2);
W2(:,[1,end])=null(D2);
[W1(2:end-1,2:end-1), L1(2:end-1)] = eig(A1, 'vector');
[W2(2:end-1,2:end-1), L2(2:end-1)] = eig(A2, 'vector');
W1([1,end],2:end-1)=G1*W1(2:end-1,2:end-1);
W2([1,end],2:end-1)=G2*W2(2:end-1,2:end-1);

[L1, L2] = ndgrid(L1, L2);
kk = sqrt(-L1-L2);

[xx,yy]=ndgrid(x);
U0=exp(-100*((xx+0.4).^2+(yy).^2)/2)+(yy.^2+xx.^2)/4;
V0=-0*(0*D*U0+0*U0*D');

UU=W1\(U0)/(W2');
VV=W1\(V0)/(W2');
VV=VV./(kk+(kk==0));
U=@(t) W1*(UU.*cos(kk*t)+VV.*sin(kk*t))*W2';

figure(1);
h=surf(xx,yy,U0);
colormap(jet(256));
camlight; shading interp; 
zlim([-1,1]); 
axis square;

dt=0.05;
nframes=1000;
for i=0:nframes-1
    set(h,'ZData',U(i*dt));
    drawnow;
end

end