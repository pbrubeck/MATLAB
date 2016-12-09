function [  ] = sectors( n )

th1=pi/3;
th2=pi/3;

a=2;
b=4;
c=1;
d=c+b-a;

[D,x]=chebD(n);

t=(x+1)/2;

z1=(a+(b-a)*t(:))*exp(1i*th1*t(:)');
z2=(c+(d-c)*t(:))*exp(1i*th2*t(:)');

z2=z2*exp(1i*(th1-th2+pi))-(d+a)*exp(1i*(th1+pi));

figure(1);
mesh(real(z1), imag(z1), 0*real(z1)); hold on;
mesh(real(z2), imag(z2), 0*real(z2)); hold off;
colormap([0,0,0]);
view(2); axis square;
end

