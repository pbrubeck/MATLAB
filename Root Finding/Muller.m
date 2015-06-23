function x2 = Muller(f, x0, x1, x2)
% Solves f(x)=0, requires three initial aproximations.
y=f([x0, x1, x2]);
while(abs(y(3))>2*eps)
    q=(x2-x1)/(x1-x0);
    q2=q*q;
    r=1+q;
    a=y*[q2, q2, 0; -q*r, -r*r, 0; q, 2*q+1, r];
    sgn=2*(a(2)>=0)-1;
    temp=x2;
    x2=x2-(x2-x1)*2*a(3)/(a(2)+sgn*sqrt(a(2)*a(2)-4*a(1)*a(3)));
    x0=x1;
    x1=temp;
    y(1)=y(2);
    y(2)=y(3);
    y(3)=f(x2);
end