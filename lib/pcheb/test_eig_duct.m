%%
xnh  = 100;
R1   = 0.3;
R2   = 1;
N    = 20;
Neig = 12;
Nrow = 3;
Ncol = ceil(Neig/Nrow);

xbox = R1*[  1 -1  -1  1  1 ];
ybox = R1*[ -1 -1   1  1 -1 ];

flag = 1;  % whether to use prescribed contour levels

lev = -0.9:0.2:0.9;  % prescribed levels, suppressing zero-contours

[lambda,x,Modes,xx,C] = solve_eig_duct(N,R1,R2,xnh,Neig);

for j = 1:Neig,
    
    subplot(Nrow,Ncol,j)
    
    Uj = expand2d(Modes(:,:,j), C);
    
    if flag
        Uj = Uj/max(max(abs(Uj)));
        contour(xx,xx,Uj,lev),  axis_sq,  xlabel x,  ylabel y
    else
        contour_sq(xx,Uj)
    end
    
    hold on
    plot(xbox,ybox,'k-')
    hold off
    
    title(['\lambda  = ',num2str(lambda(j),'%9.6f')])
    
end
