function [F] = curvedquad(z,r)
a=z([1,2,1,3]);
b=z([3,4,2,4]);
c=(b+a)/2;
d=(b-a)/2;
h=sqrt(abs(r.^2-abs(d).^2));
w=atan2(sign(r).*abs(d), h);
s=1i*sign(r).*d./abs(d);
f=cell(4,1);
for i=1:4
    if isinf(r(i))
        f{i}=@(t) c(i)-t*d(i);
    else
        f{i}=@(t) c(i)+s(i)*(abs(r(i))*exp(1i*w(i)*t)-h(i));
    end
end
F=@(x,y) (1+x)/2.*f{1}(y) + (1-x)/2.*f{2}(y) + (1+y)/2.*f{3}(x) + (1-y)/2.*f{4}(x) + ...
-1/4*( z(1)*(1+x).*(1+y) + z(2)*(1-x).*(1+y) + z(3)*(1+x).*(1-y) + z(4)*(1-x).*(1-y));
end

