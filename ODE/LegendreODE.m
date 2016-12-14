function [] = LegendreODE(n, m)
[D, x]=legD(n);
A=-diag(1-x.^2)*D*D+diag(2*x)*D;
[P,L]=eig(A, 'vector');
[L,id]=sort(L);
P=P(:,id);

% Normalization
P=bsxfun(@rdivide, P(:,m), P(end,m));
plot(acos(x),P); xlim([0,pi]);
disp(L(m));
end