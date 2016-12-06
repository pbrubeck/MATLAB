m=0;
rlim=100;
r=linspace(-rlim,rlim,2048)';
[tau,w]=gauleg(-pi,pi,256);
f=exp(1i*(m*ones(size(r))*tau-r*sin(tau)));
J=1/(2*pi)*f*w';

gold=besselj(m,r);
disp(norm(J-gold));
subplot(2,1,1); plot(r, real(J-gold));
subplot(2,1,2); plot(r, BesselJ(m,r)-gold);