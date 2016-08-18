t=linspace(0,100,2^20);
z=Zeta(0.5+1i*t);
plot3(real(z),imag(z),t);
plot(t,abs(z).^2);