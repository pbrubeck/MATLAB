function S=tubular(C, r, m)
% Generates a tubular surface from a curve using the Frenet-Serret frame.
[~,~,~,N,B]=acurve(C);
v=2*pi*(0:m-1)'/m;
S=kron(C{0}', ones(m,1))+kron(N{0}', r*cos(v))+kron(B{0}', r*sin(v));
S=reshape(S, [m, size(C{0}')]);
end