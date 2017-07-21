%DEMO_ELLIPSOIDWAVE_FIGS  plot of eave functions for first eigenmodes of a tri-axial ellipsoid
%
% This example computes first 15 eigenmodes of an ellipsoid with semi-axes 1, 1.5, and 2
% and plots the eigenvwaves. The related Helmholtz equation separated in ellipsoidal coordinates is
% discretized with the Chebyshev collocation and solved as a three-parameter eigenvalue problem
%
% See also: ELLIPSOIDWAVE_MEP, ELLIPSOID_EIGS, DEMO_ELLIPSOIDWAVE.

% Reference: M. Willatzen and L. C. Lew Yan Voon, Numerical implementation of the ellipsoidal 
% wave equation and application to ellipsoidal quantum dots, Comput. Phys. Commun. 171 (2005) 1-18.

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

fs = 26; % font size

% number of collocation nodes
n1 = 15; n2 = 15;  n3 = 15;
% semi-axes
x0 = 1;  y0 = 1.5; z0 = 2;

a2 = z0^2-x0^2; b2 = z0^2-y0^2; c = a2/b2;

neig = 5;
MX = []; MY = []; MZ = []; Mlambda = []; Mmu= []; Meta = []; Mrst = [];
RST = [0 0 0;1 0 0;0 1 0;0 0 1;1 1 0;1 0 1;0 1 1; 1 1 1];
for izb = 1:8
    fprintf('Computing first %d eigenmodes of type (r,s,t) = (%d,%d,%d) \n',neig,RST(izb,1),RST(izb,2),RST(izb,3))
    [A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,t1,t2,t3] = ellipsoidwave_mep(n1,n2,n3,x0,y0,z0,RST(izb,1),RST(izb,2),RST(izb,3));
    [lambda, mu, eta, X, Y, Z] = threepareigs(A1,B1,C1,D1,A2,B2,C2,D2,A3,B3,C3,D3,neig);
    X = [zeros(1,neig); X]; % reconstruction of b.c. in x10
    for k = 1:neig
        if real(X(end,k))<0, X(:,k) = - X(:,k); end
        if real(Y(end,k))<0, Y(:,k) = - Y(:,k); end
        if real(Z(end,k))<0, Z(:,k) = - Z(:,k); end
    end
    MX = [MX X]; MY = [MY Y]; MZ = [MZ Z];
    Mlambda = [Mlambda; lambda]; 
    Mmu = [Mmu; mu]; 
    Meta = [Meta; eta]; 
    Mrst = [Mrst; ones(neig,1)*RST(izb,:)];
end

[Meta, ord] = sort(Meta);
Mlambda = Mlambda(ord); 
Mmu = Mmu(ord);
MX = MX(:,ord); MY = MY(:,ord); MZ = MZ(:,ord);
Mrst = Mrst(ord,:);

xi1 = t1*b2;
xi2 = t2*b2;
xi3 = t3*b2;

reft1 = linspace(z0^2/b2,c,400).';
reft2 = linspace(c,1,400).';
reft3 = linspace(1,0,600).';

RefX = interp1(t1,MX,reft1,'spline');
RefY = interp1(t2,MY,reft2,'spline');
RefZ = interp1(t3,MZ,reft3,'spline');

joinx = sqrt(b2*[reft1(1:end-1);reft2;reft3(2:end)]);

MW = [];
for j = 1:length(Mlambda)
    RefZ(:,j) = RefZ(:,j)/norm(RefZ(:,j),'inf');
     if abs(RefZ(1,j))>1e-8
         alfa = RefY(end,j)/RefZ(1,j);
     else
         alfa = -RefY(end-1,j)/(reft2(end-1)-reft2(end))/RefZ(2,j)*(reft3(2)-reft3(1));
     end
     RefZ(:,j) = alfa*RefZ(:,j);
     if abs(RefX(end,j))>1e-8
         beta = RefY(1,j)/RefX(end,j);
     else
         beta = -RefY(2,j)/(reft2(2)-reft2(1))/RefX(end-1,j)*(reft1(end-1)-reft1(end));
     end
     RefX(:,j) = beta*RefX(:,j);
     
     RefX(:,j) = RefX(:,j).*abs( (reft1.^(Mrst(j,1)/2)) .* ((reft1-1).^(Mrst(j,2)/2)) .* ((reft1-c).^(Mrst(j,3)/2)) );
     RefY(:,j) = RefY(:,j).*abs( (reft2.^(Mrst(j,1)/2)) .* ((reft2-1).^(Mrst(j,2)/2)) .* ((reft2-c).^(Mrst(j,3)/2)) );
     RefZ(:,j) = RefZ(:,j).*abs( (reft3.^(Mrst(j,1)/2)) .* ((reft3-1).^(Mrst(j,2)/2)) .* ((reft3-c).^(Mrst(j,3)/2)) );
 
     if (abs(RefX(end,j))<1e-8) && (RefX(end-2,j)*RefY(2,j)>0)
          RefX(:,j) = -RefX(:,j);
     end
     if (abs(RefZ(1,j))<1e-8) && (RefZ(2,j)*RefY(end-1,j)>0)
          RefZ(:,j) = -RefZ(:,j);
     end
  
     MW(:,j) = [RefX(1:end-1,j);RefY(:,j);RefZ(2:end,j)];
     MW(:,j) = MW(:,j)/norm(MW(:,j),'inf');
     if MW(end,j)<-1e-8
         MW(:,j) = -MW(:,j);
     end
     if (abs(MW(end,j))<1e-8) && (MW(end-2,j)*MW(2,j)<0)
          MW(:,j) = -MW(:,j);
     end   
end

close all

group1 = [6 9 11 14];
group2 = [3 5 10 15];
group3 = [1 2 4 7 8 12 13];
scale1 = [0.5 1.2 0.5 1];
scale2 = [1 1 1 1];
scale3 = [1.3 1  1 0.05  1.7 0.3 1];

% Group 1
figure('Position',[100 100 1000 600]) 
plot(sqrt(b2*reft3),RefZ(:,group1)*diag(scale1),'LineWidth',3);
set(gca,'FontSize',fs)
axis([0 sqrt(b2) -2.5 2.5])
set(gca,'XTick',[0 1 sqrt(b2)])
set(gca,'XTickLabel',{'0','1','b'})

figure('Position',[100 100 450 600]) 
plot(sqrt(b2*reft2),RefY(:,group1)*diag(scale1),'LineWidth',3);
set(gca,'FontSize',fs)
axis([sqrt(b2) sqrt(a2) 0 0.6])
set(gca,'XTick',[sqrt(b2) sqrt(a2)])
set(gca,'XTickLabel',{'b','a'})
set(gca,'YTick',[0. 0.2 0.4 0.6]);

figure('Position',[100 100 500 600]) 
plot(sqrt(b2*reft1),RefX(:,group1)*diag(scale1),'LineWidth',3);
set(gca,'FontSize',fs)
axis([sqrt(a2) 2 -0.035 0])
set(gca,'XTick',[sqrt(a2) 2])
set(gca,'XTickLabel',{'a','2'})
set(gca,'YTick',[-0.03:0.01:0])
legend('6','9','11','14','Location','NorthEastOutside')

% Group 2
figure('Position',[100 100 1000 600]) 
plot(sqrt(b2*reft3),RefZ(:,group2)*diag(scale2),'LineWidth',3);
set(gca,'FontSize',fs)
axis([0 sqrt(b2) -2.2 0.5])
set(gca,'XTick',[0 1 sqrt(b2)])
set(gca,'XTickLabel',{'0','1','b'})
set(gca,'YTick',[-2 -1 0]);

figure('Position',[100 100 450 600]) 
plot(sqrt(b2*reft2),RefY(:,group2)*diag(scale2),'LineWidth',3);
set(gca,'FontSize',fs)
axis([sqrt(b2) sqrt(a2) 0 0.22])
set(gca,'XTick',[sqrt(b2) sqrt(a2)])
set(gca,'XTickLabel',{'b','a'})
set(gca,'YTick',[0 0.1 0.2])

figure('Position',[100 100 500 600]) 
plot(sqrt(b2*reft1),RefX(:,group2)*diag(scale2),'LineWidth',3);
set(gca,'FontSize',fs)
axis([sqrt(a2) 2 0 0.2])
set(gca,'XTick',[sqrt(a2) 2])
set(gca,'XTickLabel',{'a','2'})
set(gca,'YTick',[0 0.1 0.2])
legend('3','5','10','15','Location','NorthEastOutside')

% Group 3
figure('Position',[100 100 1000 600]) 
plot(sqrt(b2*reft3),RefZ(:,group3)*diag(scale3),'LineWidth',3);
set(gca,'FontSize',fs)
axis([0 sqrt(b2) -0.50 2])
set(gca,'XTick',[0 1 sqrt(b2)])
set(gca,'XTickLabel',{'0','1','b'})
set(gca,'YTick',[0 1 2]);

figure('Position',[100 100 450 600]) 
plot(sqrt(b2*reft2),RefY(:,group3)*diag(scale3),'LineWidth',3);
set(gca,'FontSize',fs)
axis([sqrt(b2) sqrt(a2) -0.2 0.8])
set(gca,'XTick',[sqrt(b2) sqrt(a2)])
set(gca,'XTickLabel',{'b','a'})
set(gca,'YTick',[0 0.4 0.8])

figure('Position',[100 100 500 600]) 
plot(sqrt(b2*reft1),RefX(:,group3)*diag(scale3),'LineWidth',3);
set(gca,'FontSize',fs)
axis([sqrt(a2) 2 -0.15 0.4])
set(gca,'XTick',[sqrt(a2) 2])
set(gca,'XTickLabel',{'a','2'})
set(gca,'YTick',[0 0.2 0.4])
legend('1','2','4','7','8','12','13','Location','NorthEastOutside')



