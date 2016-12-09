function [] = schurDecomp( N )
% Schur Decompositon Method for the Poisson equation on the L-shaped
% membrane

F1=-ones(N-2);
F2=-ones(N-2);
F3=zeros(N-2,1);

[D,x]=chebD(N);
D2=D*D;

A=D2(2:N-1,2:N-1);
[V,L]=eig(A,'vector');
[Lx,Ly]=ndgrid(L);
LL=Lx+Ly;

function u=poissonSquare(F)
    u=V*((V\F/V')./LL)*V';
end

function rhs=schurComplement(u3)
    T1=poissonSquare(u3*D2(2:N-1,1)');
    T2=poissonSquare(u3*D2(2:N-1,N)');
    rhs=(D(end)-D(1))*u3-reshape(T2*D(N,2:N-1)'-T1*D(1,2:N-1)', size(u3));
end

v1=poissonSquare(F1);
v2=poissonSquare(F2);
rhs=F3-reshape(v2*D(N,2:N-1)'-v1*D(1,2:N-1)', size(F3));

u3=cgs(@schurComplement, rhs, 1e-11, N);

u1=poissonSquare(F1-u3*D2(2:N-1,1)');
u2=poissonSquare(F2-u3*D2(2:N-1,N)');


% plot solution
uu=zeros(N,2*N-1);
uu(2:end-1,2:end-1)=[u2, u3, u1];
[xx, yy]=ndgrid(x,[x+1; x(2:end)-1]);
surf(xx, yy, uu);
zrange=max(uu(:))-min(uu(:));
xrange=max(xx(:))-min(xx(:));
yrange=max(yy(:))-min(yy(:));
daspect([1 1 2*zrange/hypot(xrange,yrange)]);
xlabel('x'); ylabel('y');
shading interp; colormap(jet(256));
end