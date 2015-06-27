function J = adaptQuad(f, a, b, x, w1, w2)
% Integrates f(x) over [a, b] using an adaptive-nested quadrature
m=(b+a)/2;
dx=(b-a)/2;
y=f(m+x*dx)*dx;
J2=y*w2(:);
J1=y(1:2:end)*w1(:); 
if(abs(J1-J2)>eps)
    J=adaptQuad(f, a, m, x, w1, w2)+adaptQuad(f, m, b, x, w1, w2);
else
    J=J2;    
end
end