%DEMO_LAME_CONV  convergence of Chebyshev collocation for Lame system
%
% In this example we compute many eigenvalues of Lame's system for
% increasing number of collocation nodes and show the convergence
% graphically
%
% If you do not have enough memory, edit file and lower neig, Nmax, and Nexact, 
% e.g., set neig=300 and Nmax=70, Nexact=100 
%
% See also: LAME_MEP, DEMO_LAME, DEMO_LAME_EXTENDED

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 14.04.2015 : fixed bug relN -> Nexact
% BP 26.08.2015 : use smaller tol and Krylov-Schur
% Last revision: 8.9.2015

neig  = 700;  % must be divided by 20
Nmax = 110;   % maximum number of nodes (we start with 30)
Nexact = 150; % relevant point (we take this as "exact" eigenvalues)

neigplus = neig + 50;
N = [30:10:Nmax Nexact]; % number of nodes to compare
 
opts = [];
opts.tol = 1e-14;
disp('Please be patient, this computation will take some time...')

values = [0.04021 0.5]; % parameter chi for the equation

for ind = 1:length(values)
    disp('------------------------------------------------------')
    fprintf('Computing eigenvalues for chi = %f for figure %d\n',values(ind),ind);
    disp('-------------------------------------------------------')
    BigMu = [];
    for k = 1:length(N)
        [A1,B1,C1,A2,B2,C2] = lame_mep(values(ind),N(k),N(k));
        A1=A1+5*B1; A2=A2+5*B2;
        fprintf('Computing %d eigenvalues for N=%d ... ',neigplus,N(k))
        tic; [lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,neigplus,opts); t1 = toc;
        fprintf('elapsed time: %f\n',t1)
        [mu,ord] = sort(mu);
        BigMu = [BigMu mu(1:neig)];
    end
    
    figure('Position',[0 0 1500 1000])
    RelDif = abs(BigMu(:,1:end-1)-BigMu(:,end)*ones(1,length(N)-1))./abs(BigMu(:,end)*ones(1,length(N)-1));
    
    RelDif20 = [];
    for j=1:neig/20
        for k=1:length(N)-1
            RelDif20(j,k) = max(RelDif(((j-1)*20 + 1):(j*20),k));
        end
    end
    
    [m1,m2] = size(RelDif20);
    
    M20 = min(zeros(m1,m2),max(-10*ones(m1,m2),log10(RelDif20+eps))); % we cut data below 1e-10 and above 1
    pcolor([M20 NaN(m1,1); NaN(1,m2+1)]);
    colorbar
    set(gca,'FontSize',30);
    set(gca,'XTick',1.5:1:length(N))
    set(gca,'XTickLabel',{'30','40','50','60','70','80','90','100','110','120'})
    set(gca,'YTick',[1 6 11 16 21 26 31 36])
    set(gca,'YTickLabel',{'0','100','200','300','400','500','600','700'})
    colorbar('YTickLabel',{'1e-10','1e-9','1e-8','1e-7','1e-6','1e-5','1e-4','1e-3','1e-2','1e-1','1e-0'},'FontSize',30)
end


