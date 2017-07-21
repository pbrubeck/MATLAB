%DEMO_WEBER_FIGS  Weber system is solved as a two-parameter eigenvalue problem 
% with Chebyshev collocation and eigenmodes are plotted
%
% This example reconstructs eigenmodes in Figure 4 (and more) (page 261) 
% from Willatzen & Lew Yan Voon,  Theory of acoustic eigenmodes in parabolic
% cylindrical enclosures, Journal of Sound and Vibration, (2005) 251--264.
%
% See also: DEMO_WEBER, WEBER_MEP, TWOPAREIGS

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 14.04.2015 : change colormap to use the same colors as in Matlab 2012
% BP 26.08.2105 : reset opts
% Last revision: 8.9.2015

N = 60;
close all % closes all open figures
mu0 = 1;
nu0 = 1;

opts = [];
opts.refine = 4; % we have to do more refine steps for accurate results

% modes with odd N
% -------------------------------------------------------------
[A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = weber_mep(N,N,1,mu0,nu0);
% we use shift to make A1 and A2 nonsingular
A1 = A1+5*B1; 
A2 = A2+5*B2;
[lambdaA,muA,XA,YA] = twopareigs(A1,B1,C1,A2,B2,C2,12,opts);
XRA = recover_bc(XA,G1,k1,r1);
YRA = recover_bc(YA,G2,k2,r2);

% N is odd, we extend v and solution YR
v = [z2; -z2(end-1:-1:1)];
YRA = [YRA; -YRA(end-1:-1:1,:)];

[tmp, ordA] = sort(abs(real(muA)));
oddmodes = [lambdaA(ordA)-5 muA(ordA)];

% mesh values for figures
matXA = 1/2*((z1.^2)*ones(1,2*N-1)-ones(N,1)*(transpose(v).^2)); 
matYA = z1*v';

% modes with even N
% -------------------------------------------------------------
[A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = weber_mep(N,N,2,mu0,nu0);
% we use shift to make A1 and A2 nonsingular
A1 = A1+5*B1; 
A2 = A2+5*B2;
[lambdaB,muB,XB,YB] = twopareigs(A1,B1,C1,A2,B2,C2,12,opts);
XRB = recover_bc(XB,G1,k1,r1);
YRB = recover_bc(YB,G2,k2,r2);

% N is even, we extend v and solution YR
v = [z2; -z2(end-1:-1:1)];
YRB = [YRB; YRB(end-1:-1:1,:)];

[tmp, ordB] = sort(abs(real(muB)));
evenmodes = [lambdaB(ordB)-5 muB(ordB)];

% mesh values for figures
matXB = 1/2*((z1.^2)*ones(1,2*N-1)-ones(N,1)*(transpose(v).^2)); 
matYB = z1*v';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot first 9 nontrivial modes in view(135,40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[0 0 1200 800]) 
oddpos = [1 3 6 7];
evenpos = [2 4 5 8 9];
fig = 1;
k = 1;
% plot odd modes 1 to 12 (with nonnegative alpha)
while fig<5
    if real(lambdaA(ordA(k))-5)>-0.1 && real(muA(ordA(k)))<-0.1
        matZ1 = real(XRA(:,ordA(k))*YRA(:,ordA(k))');
        matZ1 = matZ1 / max(max(abs(matZ1)));
        subplot(3,3,oddpos(fig));
        surf(matXA,matYA,matZ1)
        shading interp
        view(135,40)
        % title(sprintf('Odd mode %1d (alpha,beta)=(%12.8f,%12.8f)',k,lambda(ord(k))-5,mu(ord(k))))
        set(gca,'FontSize',12)
        fig = fig + 1;
    end
    k = k + 1;
end

fig = 1;
k = 1;
% plot even modes 1 to 12 (with nonnegative alpha)
while fig<6
    if real(lambdaB(ordB(k))-5)>-0.1 && real(muB(ordB(k)))<-0.1
        matZ1 = real(XRB(:,ordB(k))*YRB(:,ordB(k))');
        matZ1 = matZ1 / max(max(abs(matZ1)));
        subplot(3,3,evenpos(fig));
        surf(matXB,matYB,matZ1)
        shading interp
        view(135,40)
        % title(sprintf('Odd mode %1d (alpha,beta)=(%12.8f,%12.8f)',k,lambda(ord(k))-5,mu(ord(k))))
        set(gca,'FontSize',12)
        fig = fig + 1;
    end
    k = k + 1;
end
colormap jet  % to keep Matlab 2012 colors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot first 9 nontrivial modes in view(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',[0 0 1200 800]) 
fig = 1;
k = 1;
% plot odd modes 1 to 12 (with nonnegative alpha)
while fig<5
    if real(lambdaA(ordA(k))-5)>-0.1 && real(muA(ordA(k)))<-0.1
        matZ1 = real(XRA(:,ordA(k))*YRA(:,ordA(k))');
        matZ1 = matZ1 / max(max(abs(matZ1)));
        subplot(3,3,oddpos(fig));
        surf(matXA,matYA,matZ1)
        shading interp
        view(2)
        axis off
        % title(sprintf('Odd mode %1d (alpha,beta)=(%12.8f,%12.8f)',k,lambda(ord(k))-5,mu(ord(k))))
        set(gca,'FontSize',15)
        fig = fig + 1;
    end
    k = k + 1;
end

fig = 1;
k = 1;
% plot even modes 1 to 12 (with nonnegative alpha)
while fig<6
    if real(lambdaB(ordB(k))-5)>-0.1 && real(muB(ordB(k)))<-0.1
        matZ1 = real(XRB(:,ordB(k))*YRB(:,ordB(k))');
        matZ1 = matZ1 / max(max(abs(matZ1)));
        subplot(3,3,evenpos(fig));
        surf(matXB,matYB,matZ1)
        shading interp
        view(2)
        axis off
        % title(sprintf('Odd mode %1d (alpha,beta)=(%12.8f,%12.8f)',k,lambda(ord(k))-5,mu(ord(k))))
        set(gca,'FontSize',15)
        fig = fig + 1;
    end
    k = k + 1;
end
colormap jet  % to keep Matlab 2012 colors
