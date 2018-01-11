function [J]=confmapquad(z0,xx,yy)
% Grid
p=polygon(z0);
f=rectmap(p, 1:4);
params=parameters(f);
xmin=min(real(params.prevertex));
xmax=max(real(params.prevertex));
ymin=min(imag(params.prevertex));
ymax=max(imag(params.prevertex));
dx=xmax-xmin;
dy=ymax-ymin;

% Jacobian determinant
zz=(xmin+dx/2*(xx+1))+1i*(ymin+dy/2*(yy+1));
J=dx/2*abs(evaldiff(f,zz)).^2;
J(1:end,1:end)=0;
end