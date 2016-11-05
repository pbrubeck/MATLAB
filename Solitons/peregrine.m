function psi=peregrine(x,t)
% Peregine soliton
psi=(1-4*(1+2i*t)./(1+4*(x.^2+t.^2))).*exp(1i*t);
end