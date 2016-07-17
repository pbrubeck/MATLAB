function x2 = Muller(f, x0, x1, x2)
% Solves f(x)=0, requires three initial aproximations.
y0=f(x0);
y1=f(x1);
y2=f(x2);
tol=1E-10;
its=0;
while(max(abs(y2(:)))>tol && its<40)
    q=(x2-x1)./(x1-x0);
    q2=q.*q;
    
    a=y0.*(q2)+y1.*(-q2-q)+y2.*(q);
    c=y2.*(q+1);
    b=a+c-y1.*(q+1);
    
    sgn=2*(b>=0)-1;
    temp=x2;
    x2=x2-(x2-x1).*(2*c)./(b+sgn.*sqrt(b.*b-4*a.*c));
    x0=x1;
    x1=temp;
    y0=y1;
    y1=y2;
    y2=f(x2);
    its=its+1;
end
end