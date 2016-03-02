function [] = FieldLines(n)
t=linspace(0,2*pi,n);
r0=0.05;
ri=[r0*cos(t); r0*sin(t)];
dV=1;
V=0;
r=zeros(2,n);
for i=1:20
    k1=dV*tangent(ri);
    k2=dV*tangent(ri+k1/2);
    k3=dV*tangent(ri+k2/2);
    k4=dV*tangent(ri+k3);
    ri=ri+(k1+2*k2+2*k3+k4)/6;
    V=V+dV;
    r(:,:,i)=ri;
end

scatter(r(1,:),r(2,:));

end

function E=field(r)
E=bsxfun(@rdivide, r, dot(r,r).^(3/2));
end

function T=tangent(r)
E=field(r);
T=bsxfun(@rdivide, E, dot(E,E));
end