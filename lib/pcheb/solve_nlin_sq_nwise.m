function [dof,x,U,xx,C] = solve_nlin_sq_nwise(Nvec,xel,xnh,func,tol)
%SOLVE_NLIN_SQ_NWISE  Iteratively solve a nonlinear PDE, progressively
%                       increasing the polynomial order N as specified:
%
%                    Uxx + Uyy = func(U), Dirichlet BCs.

N0 = 0;  U0 = 0;

    for N = Nvec,
        
        [dof,x,U,xx,C,k] = mysolve(N,N0,U0,xel,xnh,func,tol);
        
        UU = expand2d(U,C);
                
        if N0==0,
            str = '';
        else
            dU = max(max(abs(UU - UU0)));
            
            str = ['  Max change from N0=',int2str(N0), ...
                   ':  du = ',num2str(dU)];
        end
        
        disp(['N = ',int2str(N),'  [',int2str(k),' iterations]',str])
        
        N0 = N;  U0 = U;  UU0 = UU;
        
    end
    
end


function [dof,x,U,xx,C,k] = mysolve(N,N0,U0,xel,xnh,func,tol)


    a0 = xel(1);  a1 = xel(end);

    [x,xx,C,~,D2,mx,Mx] = pcheb(N,xel,xnh);  n = length(x);  I = speye(n);

    [m,M,xv,yv] = mort_sq(x,mx,Mx);

    wall = (xv==a0 | xv==a1 | yv==a0 | yv==a1);
    dom  = ~wall;

    zero = zeros(size(xv));

    u = zero;
    
    if N0 > 0,
        U0 = expand2d(U0, pbary_remesh(N0,N,xel));
        u0 = U0(:);
    else
        u0 = zero; 
    end
    

    % construct Poisson operator and forcing function
    A = kron(D2,I) + kron(I,D2); 

    % mortar BCs - enforce continuity of first derivative
    A(m,:) = 0;   A = A + M;

    % solve
    dof = length(find(dom & ~m));


    du = 1.0;  k=0;


    while du > tol
        k  = k+1;
        b  = func(u0);  b(m) = 0;
        u(dom) = A(dom,dom)\b(dom);
        du = max(abs(u - u0));  u0 = u;
    end

    % convert to matrix
    U = reshape(u,n,n);


end

