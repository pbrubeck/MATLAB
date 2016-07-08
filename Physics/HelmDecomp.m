function [Phi, A1, A2, A3] = HelmDecomp(N)
% Calculates the scalar and vector potentials of a vector field
% F=-grad(Phi)+curl(A)

% Differential operators
N(1:3)=N;
[Dx,x]=chebD(N(1)); Dxx=Dx*Dx;
[Dy,y]=chebD(N(2)); Dyy=Dy*Dy;
[Dz,z]=chebD(N(3)); Dzz=Dz*Dz;
[xxx,yyy,zzz]=ndgrid(x,y,z);

% Define vector field F(x,y)=[F1, F2, F3]
r1=xxx.^2+(yyy-0.4).^2+(zzz+0.4).^2;
r2=xxx.^2+(yyy+0.4).^2+(zzz-0.4).^2;
m1=exp(-20*r1);
m2=exp(-20*r2);

F1=-m1.*(yyy-0.4)+m2.*xxx;
F2=m1.*xxx-m2.*(zzz-0.4);
F3=m1.*(zzz+0.4)+m2.*(yyy+0.4);

% Solve Poisson vector equation div(grad(P)) = -F
L={Dxx(2:end-1,2:end-1),Dyy(2:end-1,2:end-1),Dzz(2:end-1,2:end-1)};
P1=zeros(N); P2=zeros(N); P3=zeros(N);
P1(2:end-1,2:end-1,2:end-1)=kronsolve(L,-F1(2:end-1,2:end-1,2:end-1));
P2(2:end-1,2:end-1,2:end-1)=kronsolve(L,-F2(2:end-1,2:end-1,2:end-1));
P3(2:end-1,2:end-1,2:end-1)=kronsolve(L,-F3(2:end-1,2:end-1,2:end-1));

Phi=chebDiv(P1, P2, P3, Dx, Dy, Dz);
[A1, A2, A3]=chebCurl(P1, P2, P3, Dx, Dy, Dz);

figure(1);
quiver3(xxx,yyy,zzz,F1,F2,F3); axis equal;

[G1,G2,G3]=chebGrad(-Phi,Dx,Dy,Dz);
[H1,H2,H3]=chebCurl(A1,A2,A3,Dx,Dy,Dz);
[E1,E2,E3]=deal(H1+G1, H2+G2, H3+G3);

figure(2);
quiver3(xxx,yyy,zzz,E1,E2,E3); axis equal;
end

function DF=partialD(D,F,i)
sizeDF=size(F);
sizeDF(i)=size(D,1);
order=[i:ndims(F) 1:i-1];

DF=D*reshape(permute(F, order), size(D,2), []);
DF=reshape(DF, sizeDF(order));
DF=ipermute(DF, order);
end

function [F1, F2, F3]=chebGrad(P, Dx, Dy, Dz)
F1=partialD(Dx,P,1);
F2=partialD(Dy,P,2);
F3=partialD(Dz,P,3);
end

function divF=chebDiv(F1, F2, F3, Dx, Dy, Dz)
divF=partialD(Dx,F1,1)+partialD(Dy,F2,2)+partialD(Dz,F3,3);
end

function [A1,A2,A3]=chebCurl(F1, F2, F3, Dx, Dy, Dz)
A1=partialD(Dy,F3,2)-partialD(Dz,F2,3);
A2=partialD(Dz,F1,3)-partialD(Dx,F3,1);
A3=partialD(Dx,F2,1)-partialD(Dy,F1,2);
end