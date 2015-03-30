function x = MultiJacobi(f, x0)
ii=0;
err=1;
while(err>=1E-6 && ii<30)
    x=f(x0);
    err=max(abs((x-x0)./x));
    x0=x;
    ii=ii+1;
    disp(x);
end  
end