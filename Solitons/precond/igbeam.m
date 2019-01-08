function [uu] = igbeam(xi,eta,rr,p,m1,m2,q,omega,mass)
gg=exp(-omega*rr.^2/2);
cc=real(InceC(p,m1,q,1i*xi).*InceC(p,m1,q,eta)).*gg;
ss=imag(InceS(p,m2,q,1i*xi).*InceS(p,m2,q,eta)).*gg;
cc=cc/sqrt(mass(cc,cc));
ss=ss/sqrt(mass(ss,ss));
uu=(cc+1i*ss)/sqrt(2);
end