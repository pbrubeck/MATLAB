%DEMO_WEBER  Weber system is solved as a two-parameter eigenvalue problem 
% with Chebyshev collocation
%
% This example reconstructs eigenmodes in Figure 4 (and more) (page 261) 
% from Willatzen & Lew Yan Voon,  Theory of acoustic eigenmodes in parabolic
% cylindrical enclosures, Journal of Sound and Vibration, (2005) 251--264.
%
% See also: WEBER_MEP, TWOPAREIGS, DEMO_WEBER_FIGS

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
[lambda,mu,X,Y] = twopareigs(A1,B1,C1,A2,B2,C2,12,opts);
XR = recover_bc(X,G1,k1,r1);
YR = recover_bc(Y,G2,k2,r2);

% N is odd, we extend v and solution YR
v = [z2; -z2(end-1:-1:1)];
YR = [YR; -YR(end-1:-1:1,:)];

[tmp, ord] = sort(abs(real(mu)));
oddmodes = real([lambda(ord)-5 mu(ord)]);
fprintf('\nOdd modes\n');
disp('=============================');
disp('     alpha    |       beta   ');
disp('-----------------------------');
for k=1:size(oddmodes,1)
    fprintf('%13.8f | %13.8f\n',oddmodes(k,1),oddmodes(k,2))
end
fprintf('-----------------------------\n\n');


% mesh values for figures
matX = 1/2*((z1.^2)*ones(1,2*N-1)-ones(N,1)*(transpose(v).^2)); 
matY = z1*v';

% plot odd modes 1 to 12 (with nonnegative alpha)
for k=1:12
    if real(lambda(ord(k))-5)>-0.1 && real(mu(ord(k)))<-0.1
        matZ1 = XR(:,ord(k))*YR(:,ord(k))';
        figure
        surf(matX,matY,real(matZ1))
        shading interp
        title(sprintf('Odd mode %1d (alpha,beta)=(%12.8f,%12.8f)',k,lambda(ord(k))-5,mu(ord(k))))
        set(gca,'FontSize',15)
        colormap jet  % to keep Matlab 2012 colors
    end
end

% modes with even N
% -------------------------------------------------------------
[A1,B1,C1,A2,B2,C2,z1,z2,G1,k1,r1,G2,k2,r2] = weber_mep(N,N,2,mu0,nu0);
% we use shift to make A1 and A2 nonsingular
A1 = A1+5*B1; 
A2 = A2+5*B2;
[lambda,mu,X,Y] = twopareigs_ks(A1,B1,C1,A2,B2,C2,12,opts);
XR = recover_bc(X,G1,k1,r1);
YR = recover_bc(Y,G2,k2,r2);

% N is even, we extend v and solution YR
v = [z2; -z2(end-1:-1:1)];
YR = [YR; YR(end-1:-1:1,:)];

[tmp, ord] = sort(abs(real(mu)));
evenmodes = real([lambda(ord)-5 mu(ord)]);
disp('Even modes');
disp('=============================');
disp('     alpha    |       beta   ');
disp('-----------------------------');
for k=1:size(evenmodes,1)
    fprintf('%13.8f | %13.8f\n',evenmodes(k,1),evenmodes(k,2))
end
fprintf('-----------------------------\n\n');


% mesh values for figures
matX = 1/2*((z1.^2)*ones(1,2*N-1)-ones(N,1)*(transpose(v).^2)); 
matY = z1*v';

% plot even modes 1 to 12
for k=1:12
    if real(lambda(ord(k))-5)>-0.1 && real(mu(ord(k)))<-0.1
        matZ1 = XR(:,ord(k))*YR(:,ord(k))';
        figure
        surf(matX,matY,real(matZ1))
        shading interp
        title(sprintf('Even mode %1d (alpha,beta)=(%12.8f,%12.8f)',k,lambda(ord(k))-5,real(mu(ord(k)))))
        set(gca,'FontSize',15)
        colormap jet  % to keep Matlab 2012 colors
    end
end
