function [] = setlatex(fs)
% Sets latex as default text interpreter
if(nargin<1), fs=14; end
set(groot, 'DefaultAxesFontSize', fs);
set(groot, 'DefaultTextInterpreter',          'latex');
set(groot, 'DefaultLegendInterpreter',        'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
end