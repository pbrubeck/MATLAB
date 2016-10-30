function [] = sinegordon( N )
[D,x]=chebD(N); t=1+x';
Dxx=D*D;
Dtt=D*D;
Dtt(1,:)=D(end,:);

u=zeros(N);
ut0=-0*exp(-10*x.^2);
u(:,1)=ut0;
u(:,end)=exp(-10*x.^2);
BC=Dxx*u-u*Dtt';

i=0; err=1; tol=1e-10;
while err>tol
    unew=sylvester(Dxx(2:end-1,2:end-1), -Dtt(1:end-1,1:end-1)', [ut0(2:end-1), sin(u(2:end-1,2:end-1))]-BC(2:end-1,1:end-1));
    err=norm(unew-u(2:end-1,1:end-1),inf);
    u(2:end-1,1:end-1)=unew;
    i=i+1;
end
[xx,tt]=ndgrid(x,t);
surf(xx,tt,u); colormap(jet(256)); view(2);
end

