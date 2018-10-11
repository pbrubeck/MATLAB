function [E] = vnlse2(m,n,L)

% Ansatz
%spin=0; del=0; ep=pi/4; a0=2; a1=2; a2=a1;
%spin=1; del=pi/4; ep=pi/4; a0=1; a1=2.7903; a2=a1;
%spin=2; del=pi/4; ep=pi/4; a0=0.8291; a1=3.3941; a2=a1;
%spin=2; del=0; ep=pi/4; a0=0.5767; a1=3.4560; a2=a1;
%spin=4; del=pi/3; ep=pi/4; a0=0.421566597506070; a1=2.872534677296654; a2=a1;
%spin=2; del=0; ep=5*pi/16; a0=1.2722; a1=2.2057; a2=1.3310;

spin=0; del=0; ep=pi/4; a0=2.3; a1=2; a2=a1;

c=1;

% Nonlinear potential
s=0.05;
f=@(u2) -u2/s+log(1+s*u2)/s^2;
f=@(u2) -u2.^2/2;

% Linear Hamiltonian
lam=1/2;
VL=@(r) (0*r).^2;
[rr,th,jac,M,H,U,hshuff,J1,J2]=schrodpol(m,n,L,lam,VL);

% Physical domain
xx=rr.*cos(th);
yy=rr.*sin(th);
   
ii=1:m;
jj=[1:n,1];

% Ansatz
function u0=ansatz(p0,p1,p2)
u0=(p0.^((spin+1)/2)*exp(-(xx/p1).^2-(yy/p2).^2).*...
   ((cos(ep)*xx).^2+(sin(ep)*yy).^2).^(spin/2).*...
   (cos(del)*cos(spin*th)+1i*sin(del)*sin(spin*th)));
% zz=acosh((xx+1i*yy)/c);
% xi=real(zz);
% eta=imag(zz);
% u0=p0*igbeam(xi,eta,rr,3,1,3,p1*c^2,p1,M);
end


% Energy
function [E]=energy(u)
    ju=J1*u*J2';
    u2=abs(ju).^2;
    E=real(H(u,u)+jac(:)'*f(u2(:)))/2;
end

%cost=@(a0,a1,a2) energy(ansatz(a0,a1,a2));
cost=@(a0,a1,a2) energy(ansatz(a0,a1,a2));

a=[a0;a1;a2];
vars=num2cell(a);
E=cost(vars{:});
display(E);

setlatex();
figure(1);
u=ansatz(vars{:});
hp=surf(xx(ii,jj),yy(ii,jj),abs(u(ii,jj)).^2);
shading interp;
axis square;
xlim([-L,L]/sqrt(2));
ylim([-L,L]/sqrt(2));
colormap(magma(256));
colorbar();
view(2);
title(num2str(E,'$E = %f$'))
drawnow;

% Newton-Raphson for gradient
%g=agrad(cost,length(a));
%J=ahess(cost,length(a));

it=0;
e=1e-5;
y=ones(length(a),1);
da=ones(length(a),1);
tol=1e-12;
while( abs(y'*da)>tol && it<60 )
    vars=num2cell(a);
    [y,J]=fdiff(cost,a,e,vars);
    da=-J\y;
    a=a+da;
    it=it+1;
    
    u=ansatz(vars{:});
    set(hp,'ZData',abs(u(ii,jj)).^2);
    E=real(cost(vars{:}));
    title(num2str(E,'$E = %f$'))
    drawnow;
end

display(it);
display(a);

T=2*pi;
nframes=1024;
pbeam(T,nframes,u,xx,yy,jac,M,H,U,J1,J2,f);
end


function [gradf,hessf]=fdiff(f,x,e,vars)
n=length(vars);
f1=zeros(n,1);
f2=zeros(n,1);
f11=zeros(n,n);
f12=zeros(n,n);
f21=zeros(n,n);
f22=zeros(n,n);

for i=1:n
    x0=x;
    x0(i)=x0(i)+e;
    v=num2cell(x0);
    f1(i)=f(v{:});
    x0=x;
    x0(i)=x0(i)-e;
    v=num2cell(x0);
    f2(i)=f(v{:});
end
gradf=(f1-f2)/(2*e);

for i=1:n
    for j=i:n
        x0=x;
        x0(i)=x0(i)+e;
        x0(j)=x0(j)+e;
        v=num2cell(x0);
        f11(i,j)=f(v{:});
        
        x0=x;
        x0(i)=x0(i)+e;
        x0(j)=x0(j)-e;
        v=num2cell(x0);
        f12(i,j)=f(v{:});
        
        x0=x;
        x0(i)=x0(i)-e;
        x0(j)=x0(j)+e;
        v=num2cell(x0);
        f21(i,j)=f(v{:});
        
        x0=x;
        x0(i)=x0(i)-e;
        x0(j)=x0(j)-e;
        v=num2cell(x0);
        f22(i,j)=f(v{:});
        
        f11(j,i)=f11(i,j);
        f12(j,i)=f12(i,j);
        f21(j,i)=f21(i,j);
        f22(j,i)=f22(i,j);
    end
end
hessf=(f11-f12-f21+f22)/(4*e*e);
end