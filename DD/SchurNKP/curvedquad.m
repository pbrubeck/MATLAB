function [] = curvedquad(z,r)
a=z([1,2,1,3]);
b=z([3,4,2,4]);
c=(b+a)/2;
d=(b-a)/2;
h=sqrt(abs(r.^2-abs(d).^2));
w=atan2(sign(r).*abs(d), h);
v=1i*sign(r).*d./abs(d);
f=cell(4,1);
for i=1:4
    if isinf(r(i))
        f{i}=@(t) c(i)-t*d(i);
    else
        f{i}=@(t) c(i)+v(i)*(abs(r(i))*exp(1i*w(i)*t)-h(i));
    end    
end
N=16;
t=linspace(-1,1,N)';
F=@(x,y) (1-y)/2.*f{4}(x) + (1+y)/2.*f{3}(x) + (1-x)/2.*f{2}(y) + (1+x)/2.*f{1}(y) + ...
-1/4*( z(1)*(1+x).*(1+y) + z(2)*(1-x).*(1+y) + z(3)*(1+x).*(1-y) + z(4)*(1-x).*(1-y));

[xx,yy]=ndgrid(t);


figure(1);
plot(f{1}(t), 'linewidth', 2); hold on;
plot(f{2}(t), 'linewidth', 2); 
plot(f{3}(t), 'linewidth', 2);
plot(f{4}(t), 'linewidth', 2); 
plot(F(xx,yy));
plot(F(xx,yy).');
hold off;
axis equal;

end

