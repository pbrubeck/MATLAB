function y=Gamma(z)
% Computes the gamma function of a complex number.
b=real(z)<0.5;
z=b+(1-2*b).*z-1;
g=7;
k=10;
c=zeros(2*k-1,1);
fact=2;
for j=0:k-1
    c(2*j+1)=fact*exp(j+g+0.5)*(j+g+0.5)^(-j-0.5);
    fact=fact*(2*j+1)/2;
end
p=Chebyshev(2*k-1)*c;

r=ones(size(z));
s=p(1)/2*r;
for i=2:k
    r=r.*(z-i+2)./(z+i-1);
    s=s+p(2*i-1)*r;
end
t=z+g+0.5;
y=t.^(z+0.5).*exp(-t).*s;
y(b)=pi./(y(b).*sin(-pi*z(b)));
end
