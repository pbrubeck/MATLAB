function [] = HelmElliptical(a, b, N, k)
f=sqrt(a^2-b^2);
[A1,B1,A2,B2,u,v]=chebLapEll(a,b,2*N-1,N);
A1=A1(2:end-1,2:end-1);
B1=B1(2:end-1,2:end-1);
C1=eye(2*N-3);
C2=-eye(N);

a=k^2;
lam=-2*a;
rmf=cos(pi/2*u(2:end-1)/u(end));
amf=cos(k*v(:));
[lam,a,X1,X2] = newton_mep(A1,B1,C1,A2,B2,C2,rmf,amf,lam,a);
%[lam,a,X1,X2] = twopareigs(A1,B1,C1,A2,B2,C2,k);
q=-lam*f^2/4;

rmf=zeros(N,1);
rmf(2:end,:)=bsxfun(@times, conj(X1(N-1,:)), X1(1:N-1,:));
amf=bsxfun(@times, conj(X2(1,:)), X2);
rmf=rmf/max(rmf);
amf=amf/max(amf);
u=u(1:N);

display(a);
display(q);

figure(1); plot(u,rmf);
figure(2); plot(v,amf);

xx=f*cosh(u)*cos(v);
yy=f*sinh(u)*sin(v);
zz=rmf*amf.';

figure(3);
surfl(xx(:,[1:end 1]),yy(:,[1:end 1]),zz(:,[1:end 1]),'light'); 
shading interp;
axis off;
colormap(jet(256));
zrange=max(zz(:))-min(zz(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
view(0,90);
end