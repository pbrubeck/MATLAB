function Delta = extract_regular_part_np(Delta,opts)

%EXTRACT_REGULAR_PART_NP  Regular projections of Delta matrices from singular MEP 
%
% D = EXTRACT_REGULAR_PART_NP(Delta,opts)
% returns "Delta" matrices of the projected regular problem
% with D{i} = P*Delta{i}*Q, i = 1, ..., k+1 such that D{1} is of full rank, 
% matrices (D{1})^(-1)*D{2}, ..., (D{1})^(-1)*D{k+2} commute (hopefully)
% and P and Q are maximal possible matrices with orthogonal columns. It
% uses staircase algorithm with column-row and row-column compressions.
%
% opts:
%  - oneway : 0 (set to 1 if it is enough to do just the column row reduction cycle)

% Reference: A. Muhic, B. Plestenjak: On the quadratic two-parameter 
% eigenvalue problem and its linearization, Linear Algebra Appl. 432 (2010) 2529-2542.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<2, opts=[]; end
if isfield(opts,'oneway'), oneway = opts.oneway; else oneway = 0; end

full = 0;
while ~full
    [Delta, full] = staircase_step_cr_np(Delta, opts);
end;
full = oneway;
while ~full
    [Delta, full] = staircase_step_rc_np(Delta, opts);
end;
