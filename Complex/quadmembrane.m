function [] = quadmembrane(a,b,h,N,k)
% Dirichlet eigenmodes of the laplacian on the trapezoid
vertex=[a+1i*h; -a+1i*h; -b-1i*h; b-1i*h]/2;
corners=[1 2 3 4];
rectlap(vertex, corners, N, k);
end

