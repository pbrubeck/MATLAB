%DEMO_BESSELWAVE1  Bessel wave system is solved as a two-parameter eigenvalue problem 
% with Chebyshev collocation
% 
% This example reconstructs Tables 1, 2 and 3 in  
% Lew Yan Voon & Willatzen, Helmholtz equation in parabolic rotational 
% coordinates: application to wave problems in quantum mechanics and
% acoustics, Math. Comp. Sim. 65 (2004) 337--349.
%
% See also: DEMO_BESSELWAVE2, DEMO_BESSELWAVE3, DEMO_BESSELWAVE_FIGS, BESSELWAVE_MEP, TWOPAREIGS

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

N = 60;
XiSq = 1;
EtaSq = 1;

% We compute enough modes (8 for p=0,..,8) and then select the smallest sqrt(mu)
% We select only modes with lambda>0 (since some are close to zero and
% negative, we take lambda>-1)

% -------------------------------------------------------------------
% Table 1 - Dirichlet & Dirichlet
modes = [];
for p=0:8
    [A1,B1,C1,A2,B2,C2,z1,z2,G1,G2] = besselwave_mep(N,N,p,sqrt(XiSq),sqrt(EtaSq),[1 1]);
    [lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,8);
    modes = [modes; p*ones(length(lambda),1) lambda sqrt(mu)];
end

posmodes = modes(find(real(modes(:,2))>-1),:);
[tilda,ord] = sort(real(posmodes(:,3)));

Table1 = posmodes(ord(1:20),:);  % first 20 modes
fprintf('\nFirst 20 modes for Dirichlet & Dirichlet b.c.\n');
disp('=============================================');
disp('p  |    lambda   |    omega  ');
disp('------------------------------');
for k=1:20
    fprintf('%1i  | %11.8f | %11.8f \n',posmodes(ord(k),1),posmodes(ord(k),2),posmodes(ord(k),3))
end
fprintf('------------------------------\n\n');

% -------------------------------------------------------------------
% Table 2 - Neumann & Neumann
modes = [];
for p=0:8
    [A1,B1,C1,A2,B2,C2,z1,z2,G1,G2] = besselwave_mep(N,N,p,sqrt(XiSq),sqrt(EtaSq),[2 2]);
    % we have eigenvalues mu = 0, therefore we shift mu -> mu - 10 to make
    % Delta2 and A1 and A2 nonsingular
    A1=A1+10*C1; 
    A2=A2+10*C2;
    [lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,8);
    modes = [modes; p*ones(length(lambda),1) lambda sqrt(mu-10)];
end

posmodes = modes(find(real(modes(:,2))>-1),:);
[tilda,ord] = sort(real(posmodes(:,3)));

Table2 = posmodes(ord(1:20),:);  % first 20 modes
disp('First 20 modes for Neumann & Neumann b.c.');
disp('=============================================');
disp('p  |    lambda   |    omega  ');
disp('------------------------------');
for k=1:20
    fprintf('%1i  | %11.8f | %11.8f \n',posmodes(ord(k),1),posmodes(ord(k),2),posmodes(ord(k),3))
end
fprintf('------------------------------\n\n');


% -------------------------------------------------------------------
% Table 3 - Dirichlet & Neumann
modes = [];
for p=0:8
    [A1,B1,C1,A2,B2,C2,z1,z2,G1,G2] = besselwave_mep(N,N,p,sqrt(XiSq),sqrt(EtaSq),[1 2]);
    % we have singular A1 and A2, therefore we shift lambda -> lambda - 10 to make
    % A1 and A2 nonsingular
    A1=A1+10*B1; 
    A2=A2+10*B2;
    [lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,8);
    modes = [modes; p*ones(length(lambda),1) lambda-10 sqrt(mu)];
end

posmodes = modes(find(real(modes(:,2))>-1),:);
[tilda,ord] = sort(real(posmodes(:,3)));

Table3 = posmodes(ord(1:20),:);  % first 20 modes
disp('First 20 modes for Dirichlet & Neumann b.c.');
disp('=============================================');
disp('p  |    lambda   |    omega  ');
disp('------------------------------');
for k=1:20
    fprintf('%1i  | %11.8f | %11.8f \n',posmodes(ord(k),1),posmodes(ord(k),2),posmodes(ord(k),3))
end
fprintf('------------------------------\n');
