function [max]=superposition(A, w, phi)
global w2; 
global w3;
w2=w.*w;
w3=w2.*w;
t=0;
ii=0;
while(~isnan(t))
t=maximize(A, w, phi, t);
max=real(A*transpose(exp(1i*(w*t+phi))));
tmax=t;
t=solve(A, w, phi, max+1E-10, t);
fprintf('i=%d \t max=%f \t tmax=%f \t sol=%f \n', ii, max, tmax, t);
ii=ii+1;
end
end

function t=solve(A, w, phi, max, t)
global w2; 
y=0;
ii=0;
while(ii<40 && (ii==0 || abs(y)>1E-15))
    psi=(A.*exp(1i*(w*t+phi))).';
    y=real(sum(psi))-max;
    y1=-imag(w*psi);
    y2=-real(w2*psi);
    r=y/y1;
    t=t-r/(1-r*y2/(2*y1));
    ii=ii+1;
end
if(ii==40)
    t=NaN;
end
end

function t=maximize(A, w, phi, t)
global w2; 
global w3;
y1=0;
ii=0;
while(ii<20 && (ii==0 || abs(y1)>1E-15))
    psi=(A.*exp(1i*(w*t+phi))).';
    y1=-imag(w*psi);
    y2=-real(w2*psi);
    y3=imag(w3*psi);
    r=y1/y2;
    t=t-r/(1-r*y3/(2*y2));
    ii=ii+1;
end
end