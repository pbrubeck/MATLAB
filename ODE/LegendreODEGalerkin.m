function [] = LegendreODEGalerkin(n,m)
[D,x]=legD(n);
VL=VandermondeGeg(1/2, x);
VG=VandermondeGeg(3/2, x);
M=inv(VL*VL');
S=VG\D;
K=S'*S;

[P,L]=eig(K,M, 'vector');
[L,id]=sort(L);
P=P(:,id);

% Normalization
P=bsxfun(@rdivide, P(:,m), P(end,m));
plot(acos(x),P); xlim([0,pi]);
disp(L(m));
end