%DEMO_LAME  Lame system is solved as a two-parameter eigenvalue problem 
% with Chebyshev collocation
%
% This example reconstructs Table 1 and Table 2 (page 246) from 
% Morrison & Lewis, Charge singularity at the corner of a flat plate,
% SIAM J. Appl. Math. 31 (1976) 233--250.
%
% See also: LAME_MEP, TWOPAREIGS, DEMO_LAME_EXTENDED, DEMO_LAME_CONV

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% BP 26.08.2015 : use the main routine twopareigs 
% Last revision: 8.9.2015

N = 60;

t = [1:7 9:15]/8; % values [0.125 ... 1.875] xi/pi in Table 1
value = [];

for k=1:length(t)
    [A1,B1,C1,A2,B2,C2] = lame_mep(t(k)*pi,N,N); % we discretize system with Chebyhev spectral collocation
    A1=A1+5*B1; A2=A2+5*B2;                            % we shift because A2 is singular
    [lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,1);     % we need just the smallest eigenvalue, KrylovSchur is faster from eigs for just one eigenvalue
    value(k,:) = [max(roots([1 1 -mu])) mu];                  % positive root of z^2+z-mu = 0
end 

Table1 = [t(:) value];
disp('Table 1 in Morrison and Lewis, SIAM J. Appl. Math. 31 (1976) 233--250');
disp('=====================================================================');
disp('hi/pi   | rho        | lambda');
disp('---------------------------------');
for k=1:length(t)
    fprintf('%7.5f | %10.8f | %10.8f\n',t(k),value(k,:))
end
fprintf('---------------------------------\n\n');

t = [0.04021 0.11610 0.16670 0.22432 0.28858 0.875 0.900 0.925 0.950 1.050 1.075 1.100 1.125 1.875 1.900 1.925 1.950]; % values xi/pi in Table 2
value = [];

for k=1:length(t)
    [A1,B1,C1,A2,B2,C2] = lame_mep(t(k)*pi,N,N);
    A1=A1+5*B1; A2=A2+5*B2; % we shift because A2 is singular
    [lambda,mu] = twopareigs(A1,B1,C1,A2,B2,C2,1);
    value(k,:) = [max(roots([1 1 -mu])) mu];                  % positive root of z^2+z-mu = 0
end 

Table2 = [t(:) value];
disp('Table 2 in Morrison and Lewis, SIAM J. Appl. Math. 31 (1976) 233--250');
disp('=====================================================================');
disp('hi/pi   | rho        | lambda');
disp('---------------------------------');
for k=1:length(t)
    fprintf('%7.5f | %10.8f | %10.8f\n',t(k),value(k,:))
end
fprintf('---------------------------------\n');
