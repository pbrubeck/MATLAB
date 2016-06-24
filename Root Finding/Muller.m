function x2 = Muller(f, x0, x1, x2)
% Solves f(x)=0, requires three initial aproximations.
y1=f(x0);
y2=f(x1);
y3=f(x2);
tol=1E-10;
its=0;
while(max(abs(y3(:)))>tol)
    q=(x2-x1)./(x1-x0);
    q2=q.*q;
    r=1+q;
    
    a1=y1.*(q2)+y2.*(-q.*r)+y3.*(q);
    a2=y1.*(q2)+y2.*(-r.*r)+y3.*(2*q+1);
    a3=y3.*(r);
    
    sgn=2*(a2>=0)-1;
    temp=x2;
    x2=x2-(x2-x1).*(2*a3)./(a2+sgn.*sqrt(a2.*a2-4*a1.*a3));
    x0=x1;
    x1=temp;
    y1=y2;
    y2=y3;
    y3=f(x2);
    its=its+1;
end
disp(its);
end