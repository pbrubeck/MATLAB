function [] = Lmembrane(N, k)
% Dirichlet eigenmodes of the laplacian on the L-shaped membrane
vertex=[1i, -1+1i, -1-1i, 1-1i, 1, 0];
corners=[5 1 2 4];
rectlap(vertex, corners, N, k);
end