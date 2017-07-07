function [] = gordonhall( n )
a=0.5;
b=0.5;
th=atan2(a,b);

p=[1+1i; 1i; 1i*b; a*cos(th)+1i*b*sin(th)];
q=[1+1i; a*cos(th)+1i*b*sin(th); a; 1];

f1=mapline(p(2),p(1));
f2=mapline(p(3),p(2));
f3=mapellipse(0,p(3),p(4));
f4=mapline(p(4),p(1));

g1=mapline(q(2),q(1));
g2=mapellipse(0,q(3),q(2));
g3=mapline(q(3),q(4));
g4=mapline(q(4),q(1));


F=@(x,y) (1-y)/2.*f3(x)+(1+y)/2.*f1(x)+...
    (1-x)/2.*(f2(y)-(1+y)/2.*f2(1)-(1-y)/2.*f2(-1))+...
    (1+x)/2.*(f4(y)-(1+y)/2.*f4(1)-(1-y)/2.*f4(-1));

G=@(x,y) (1-y)/2.*g3(x)+(1+y)/2.*g1(x)+...
    (1-x)/2.*(g2(y)-(1+y)/2.*g2(1)-(1-y)/2.*g2(-1))+...
    (1+x)/2.*(g4(y)-(1+y)/2.*g4(1)-(1-y)/2.*g4(-1));


[xx,yy]=ndgrid(linspace(0,1,n),linspace(-1,1,n));
figure(1); clf;
plot(F(xx,yy)); hold on;
plot(G(xx,yy).'); hold off;
axis equal;

end

