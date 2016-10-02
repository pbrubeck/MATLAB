function [] = LaplaceMap(N)
[Dx,x]=chebD(N(1)); Dxx=Dx*Dx;
[Dy,y]=chebD(N(1)); y=y(:)'; Dyy=Dy*Dy;

u1=[1+1i+1i*y; -1+1i*y];
u2=[3i/2+(1+1i/2)*x, -1i/2+(1+1i/2)*x];

RHS=-Dxx(:,[1,end])*u1-u2*Dyy(:,[1,end])';
ww=zeros(N);
ww([1,end],:)=u1;
ww(:,[1,end])=u2;
ww(2:end-1,2:end-1)=sylvester(Dxx(2:end-1,2:end-1), Dyy(2:end-1,2:end-1)', RHS(2:end-1,2:end-1));

plot(ww); hold on; plot(ww.'); hold off; axis equal;
end

