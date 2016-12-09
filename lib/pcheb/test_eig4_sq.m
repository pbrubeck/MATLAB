%% Program 39 of Trefethen (Chapter 14)
xnh  = 50;
xel  = [-1,0,1];
N    = 20;
Neig = 16;
Nrow = 4;
Ncol = ceil(Neig/Nrow);

[lambda,x,Modes,xx,C] = solve_eig4_sq(N,xel,xnh,Neig);

for j = 1:Neig,
    
    subplot(Nrow,Ncol,j)
    
    Uj = expand2d(Modes(:,:,j), C);
    
    mu = lambda(j);
    
    contour_sq(xx,Uj)
    
    title(num2str(mu,'%7.4f'))
    
end
