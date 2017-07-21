%DEMO_HYDROGEN_MP   first four modes for hydrogen molecular ion H2+ in 2D 
% for 0.1<R<6 with multiprecision refinement
% 
% This example reconstructs Table 1 and Table 2 (page 2201) 
% from Patil,  Hydrogen molecular ion and molecule in two dimensions, 
% Journal of Chemical Physics 118, (2003) 2197--2205.
% 
% See also: HYDROGEN_MODES

% Reference: B. Plestenjak, C. I. Gheorghiu, M. E. Hochstenbach: Spectral
% collocation for multiparameter eigenvalue problems arising from separable
% boundary value problems, J. Comp. Phys. 298 (2015) 585-601

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 2.1.2017

is_numeric_type_supported('mp'); % check if MCT is installed

opts.fp_type='mp';

valuesR = mp([0.1:0.1:0.5 0.51 0.515 0.52 0.53 0.6 0.7 0.8 1 1.5 2 3 4 4.5 5 6]);

disp('      |             1sg+                    |             2pu+                    |             2sg+                    |            2pu-');
disp('R     |     S                  E            |     S                  E            |     S                  E            |     S                  E            ');
disp('--------------------------------------------------------------------------------------------------------------------------------------------------------------');
for R=valuesR
    if R<=2
        b = 1;
    else
        b = 2;
    end
    [S, E] = hydrogen_modes(R,100,100,b,opts);
    fprintf('%5.3f | %16.14f  %16.14f | %16.14f  %16.14f | %16.14f  %16.14f | %16.14f  %16.14f\n',...
        [R S(1) E(1) S(2) E(2) S(3) E(3) S(4) E(4)])
end
