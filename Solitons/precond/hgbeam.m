function [uu] = hgbeam(xx,yy,m,n,omega)
c1=[zeros(1,m),1];
c2=[zeros(1,n),1];
uu=sqrt(omega)*HermitePsi(c1,sqrt(omega)*xx).*HermitePsi(c2,sqrt(omega)*yy);
end

