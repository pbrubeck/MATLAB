function [D0, D1, D2, full] = staircase_step_rc(D0, D1, D2, opts)

%STAIRCASE_STEP_RC  Row-Column reduction step of staircase algorithm for Delta matrices
%
% [D0, D1, D2, full] = STAIRCASE_STEP_RC(D0, D1, D2, opts)
% reduces D0 to full column rank using row and column compression on (D0, D1), (D0, D2)
%
% See also: STAIRCASE_STEP_CR, STAIRCASE_STEP_RC_NP, STAIRCASE_STEP_CR_NP

% Reference: A. Muhic, B. Plestenjak: On the quadratic two-parameter 
% eigenvalue problem and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% This file is no longer used in MultiParEig 2.0, but is kept for backward compatibility

% Last revision: 8.9.2015

if nargin<4, opts=[]; end

opts.call = 'RedRowCol 1';
opts.firstlastgap = 2;
[ran1, U1, tilda1, tilda2] = extended_svd(D0,opts);

% row nullity
if ran1 == size(D0,1)
    full = 1;
    return;
else
    full = 0;
end

D1 = U1'*D1;
D2 = U1'*D2;
POM1 = D1(ran1+1:end,:);
POM2 = D2(ran1+1:end,:);

opts.call = 'RedRowCol 2';
opts.firstlastgap = 1;
[ran2, tilda1, tilda2, V2] = extended_svd([POM1;POM2], opts);

if ran2 == size(POM1,2)
    % reduction leads to empty matrices
    D0 = []; D1 = []; D2 = [];
    return
end

D0 = U1(:,1:ran1)'*D0*V2(:,ran2+1:end);
D1 = D1(1:ran1,:)*V2(:,ran2+1:end);
D2 = D2(1:ran1,:)*V2(:,ran2+1:end);
