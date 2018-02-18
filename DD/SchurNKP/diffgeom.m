function [jac,G11,G12,G22] = diffgeom(F,u,v)
% Initialize audi grid
m=length(u);
n=length(v);
[u,v] = ndgrid(u,v);
[u,v] = ainit(u,v,2);

% Surface parametrization
S = [real(F(u,v)); imag(F(u,v))];

% Differential geometry
J = ajac(S);    % Jacobian
G = J'*J;       % Metric

detG=det(G);
gg=G{0};

G11=reshape(gg(1,1,:), m,n);
G12=reshape(gg(1,2,:), m,n);
G22=reshape(gg(2,2,:), m,n);
jac=reshape(sqrt(max(detG{0},0)), m,n);
end