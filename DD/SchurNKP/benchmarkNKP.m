function [] = benchmarkNKP(ms)

L=3;
R1=1;  R2=1;
z1=L;  z2=-L;
m1=1;  m2=1;
function F=binrhs(z,rho)
    r1=hypot(rho,z-z1);
    r2=hypot(rho,z-z2);
    F1=-m1*(r1<=R1).*sinc(pi*r1/R1);
    F2=-m2*(r2<=R2).*sinc(pi*r2/R2);
    F=F1+F2;
end
metric = @(z,rho) abs(rho);
RHS=@binrhs;

[z0,quads,curv]=bingrid(R1,R2,L);

% Mesh refinement
ref=1;
for j=1:ref
    [~, adj, bnd] = meshtopo(quads);
    [z0,quads,curv]=quadmeshrefine(z0,quads,curv,adj,bnd);
end

ngrid=zeros(size(ms));
pcalls=zeros(size(ms));
for i=1:numel(ms)
    m=ms(i);
    [lam, uu, X, Y, relres, calls] = meshSchurNKP(m, z0, quads, curv, 0, RHS, metric);
    delete(findall(gcf,'Type','light'));
    drawnow;
    pcalls(i)=calls;
    ngrid(i)=numel(uu);
    display(relres);
    display(calls);
end

figure(20);
semilogx(ngrid,pcalls,'--*b','Linewidth',1);
set(gcf,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('Total gridpoints');
title('Preconditioner calls by GMRES');
ylim([0,100]);
xlim([1E5,5E6]);
set(gca,'XTick',ngrid);
set(gca,'XTickLabel',num2str(ngrid(:),'%.2G'));
set(gca,'XMinorTick','off');

dp=0.15*min(pcalls);
for i=1:numel(ms)
    text(ngrid(i),pcalls(i)-dp,sprintf('$n=%d$',ms(i)),...
        'FontSize',14,'HorizontalAlignment','center');
end

end

