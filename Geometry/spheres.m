function []=spheres(n,m)
x=zeros(n,m);
y=zeros(n,m);
z=zeros(n,m);

% Naive
for i=0:n-1
    u=i/(n-1);
    theta=2*pi*u;
    for j=0:m-1
        v=j/(m-1);
        phi=pi*v;
        t=sin(phi);
        x(i+1,j+1)=t*cos(theta);
        y(i+1,j+1)=t*sin(theta);
        z(i+1,j+1)=cos(phi);
    end
end
subplot(2,2,1);
surf(x,y,z);
axis equal

x=zeros(n,m);
y=zeros(n,m);
z=zeros(n,m);

% Equal area subdivision
for i=0:n-1
    u=i/(n-1);
    theta=2*pi*u;
    for j=0:m-1
        v=j/(m-1);
        phi=acos(2*v-1);
        t=2*sqrt(v*(1-v));
        x(i+1,j+1)=t*cos(theta);
        y(i+1,j+1)=t*sin(theta);
        z(i+1,j+1)=2*v-1;
    end
end
subplot(2,2,2);
surf(x,y,z);
axis equal

x=zeros(n,m);
y=zeros(n,m);
z=zeros(n,m);

% Cartesian grid mapping
for i=0:n-1
    u=2*(mod(i,n/2))/(n/2-1)-1;
    s=2*(2*i>=n)-1;
    for j=0:m-1
        v=2*j/(m-1)-1;
        if(u==0 && v==0)
            r=0; theta=0;
        elseif(u*u>v*v)
            r=u; theta=pi/4*(v/u);
        else
            r=v; theta=pi/4*(2-u/v);
        end
        phi=acos(1-r^2);
        t=r*sqrt(2-r^2);
        x(i+1,j+1)=s*t*cos(theta);
        y(i+1,j+1)=t*sin(theta);
        z(i+1,j+1)=s*(1-r^2);
    end
end
subplot(2,2,3);
surf(x,y,z);
axis equal
end