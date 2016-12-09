%% Program 23a of Trefethen (Chapter 9) - unperturbed Laplacian
xnh  = 50;
xel  = [-1,0,1];
N    = 12;
Neig = 6;
Nrow = 2;
Ncol = ceil(Neig/Nrow);

[lambda,x,Modes,xx,C] = solve_eig_sq(N,xel,xnh,Neig);

for j = 1:Neig,
    
    subplot(Nrow,Ncol,j)
    
    Uj = expand2d(Modes(:,:,j), C);
    
    mu = (4/pi^2)*lambda(j);
    
    contour_sq(xx,Uj)
    
    title(['(4/\pi^2)\lambda  = ',num2str(mu,'%12.9f')])
    
end
