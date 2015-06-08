function [] = lollipop(w)
% Displays lollipop plot of complex signal
stem([real(w(:)) imag(w(:))]);
end

