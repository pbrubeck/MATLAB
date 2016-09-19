function [Delta0,Delta1,Delta2] = extract_regular_part(Delta0,Delta1,Delta2,opts)

%EXTRACT_REGULAR_PART  Regular projections of Delta matrices from singular 2EP 
%
% [D0,D1,D2] = EXTRACT_REGULAR_PART(Delta0,Delta1,Delta2,tol)
% returns "Delta" matrices D0,D1,D2 of the projected regular problem
% with D0 = P*Delta0*Q, D1 = P*Delta1*Q, D2 = P*Delta2*Q such that 
% matrices (D0)^(-1)*D1 and (D0)^(-1)*D2 commute and P and Q are maximal 
% possible matrices with orthogonal columns.
% Matrices Delta arise from a two-parameter eigenvalue problem
%
% opts:
%  - oneway : 0 (set to 1 if it is enough to do just the column row reduction cycle)
%
% See also: EXTRACT_REGULAR_PART_NP

% Reference: A. Muhic, B. Plestenjak: On the quadratic two-parameter 
% eigenvalue problem and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% This file is no longer used in MultiParEig 2.0, but is kept for backward compatibility

% Last revision: 8.9.2015

if nargin<4, opts=[]; end
if isfield(opts,'oneway'), oneway = opts.oneway; else oneway = 0; end

full = 0;
while ~full
    [Delta0, Delta1, Delta2, full] = staircase_step_cr(Delta0, Delta1, Delta2, opts);
end;
full = oneway;
while ~full
    [Delta0, Delta1, Delta2, full] = staircase_step_rc(Delta0, Delta1, Delta2, opts);
end;
