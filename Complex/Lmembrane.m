function [lam] = Lmembrane(N, k)
% Dirichlet eigenmodes of the laplacian on the L-shaped membrane
vertex=[1; 1+1i; -1+1i; -1-1i; -1i; 0];
corners=[1; 2; 4; 5];
lam=rectlap(vertex, corners, N, k);
end