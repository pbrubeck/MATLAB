function out = linsystem(A, x0, in1, in2)
[P,D]=eig(A);
d=diag(D);
PC=P*diag(P\x0);
if nargin==3
    out=PC*exp(d*in1);
else
    t=2;
    it=0;
    err=1;
    PCD=PC*D;
    B=[PC(in2,:); PCD(in2,:)];
    while(abs(err)>1E-15 && it<20)
        f=B*exp(d*t);
        err=(f(1)-in1)/f(2);
        t=t-err;
        it=it+1;
    end
    out=t;
end

t=linspace(0, 50, 100);
xt=PC*exp(d*t);
plot(t,xt);

end