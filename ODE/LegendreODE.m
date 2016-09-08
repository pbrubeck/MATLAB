function [] = LegendreODE(n, m)
[D, x, w]=legD(n);
OP=-diag(1-x.^2)*D*D+diag(2*x)*D;
[P,L]=eig(OP, 'vector');
[L,id]=sort(L);
P=P(:,id);

% Normalization
P=normc(P, w);
P=bsxfun(@times, P(:,m), sign(P(end,m)).*sqrt(2./(2*m-1)));
plot(acos(x),P); xlim([0,pi]);
disp(L(m));
end