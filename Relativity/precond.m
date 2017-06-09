function [uu,res,its] = precond(eqn, green, rhs, u0, maxit, tol)
% bicgstab helper for 2D elliptic equations
afun=@(x) reshape(eqn(reshape(x, size(u0))), [], 1);
pfun=@(x) reshape(green(reshape(x, size(u0))), [], 1);
[x,~,res,its]=bicgstab(afun,rhs(:),tol,maxit,pfun,[],u0(:));
%[x,~,res,its]=gmres(afun,rhs(:),[],tol,maxit,pfun,[],u0(:));
uu=reshape(x,size(u0));
end