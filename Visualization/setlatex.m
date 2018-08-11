function [] = setlatex()
% Sets latex as default text interpreter
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultTextInterpreter',          'latex');
set(groot, 'DefaultLegendInterpreter',        'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
end