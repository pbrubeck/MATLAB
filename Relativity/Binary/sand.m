n=128;
x=linspace(0,1,n);
y=linspace(0,2,n);
[xx,yy]=meshgrid(x,y);
d=sqrt(xx);
ww=(2 - 6*(yy./d).^2 + 4*(yy./d).^3)./d;
ww(yy>d)=0;
figure(1);
contourc(x,y,ww,16);


X=[-37; -40; -45; -50; -55; -67];
Y=[-68; -74.1; -81; -89; -97.4; -117];
A=[ones(size(X)), X];
Z=(A'*A)\A'*Y;
plot(X,Y,X,A*Z);

Td=[-105; -50; 0; 25; 120; 150];
To=[ones(size(Td)), Td]*Z;

display([Td, To])
