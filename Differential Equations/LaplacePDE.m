function u = LaplacePDE(g, L, n)
%LAPLACEPDE Summary of this function goes here
%   Detailed explanation goes here
k=1:128;
x=linspace(0,L,n);
y=linspace(0,L,n);
X=sin(pi/L*k(:)*x);
Y=sinh(pi/L*k(:)*y);
c=sineSeries(g, L, 1024);
c=c(k);
r=csch(pi*k(:));
u=zeros(n,n);
for i=1:n
    for j=1:n
        u(n-j+1,i)=c*(X(:,i).*Y(:,j).*r);
    end
end
colormap(jet(64));
image(u,'CDataMapping','scaled');
colorbar();
end