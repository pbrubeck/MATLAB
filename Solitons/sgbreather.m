function [u, v] = sgbreather(w,x,t)
% sine-Gordon breather
c=sqrt(1-w^2);
u=4.*atan(c/w*cos(t*w).*sech(c*x));
v=4*c*w*sech(c*x)./(w^2+c^2*(cos(t*w).*sech(c*x)).^2).* ...
  (c*cos(t*w).*tanh(c*x)-w*sin(t*w));
end

