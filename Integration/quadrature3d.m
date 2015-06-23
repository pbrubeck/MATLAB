function J = quadrature3d(f, x, y, z, t, w)
% Integrates f over the 3D region bounded by the functions x, y(x), z(x,y)
% using a standard quadrature rule for the [-1, 1] interval.
n=length(t);
n3=n*n*n;
ndx=zeros(1, n3);
ndy=zeros(1, n3);
ndz=zeros(1, n3);
wgt=zeros(1, n3);
x0=(x(2)+x(1))/2;
dx=(x(2)-x(1))/2;
for i=1:n
    xi=x0+dx*t(i);
    lim=y(xi); a=lim(1); b=lim(2);
    y0=(b+a)/2;
    dy=(b-a)/2;
    dA=w(i)*dx*dy;
    ndx((i-1)*n*n+1:i*n*n)=xi;
    for j=1:n
        yj=y0+dy*t(j);
        lim=z(xi, yj); a=lim(1); b=lim(2);
        z0=(b+a)/2;
        dz=(b-a)/2;
        dV=w(j)*dA*dz;
        start=(i*n-n+j-1)*n+1;
        finish=(i*n-n+j)*n;
        ndy(start:finish)=yj;
        ndz(start:finish)=z0+dz*t;
        wgt(start:finish)=w*dV;
    end
end
J=f(ndx, ndy, ndz)*wgt(:);
end