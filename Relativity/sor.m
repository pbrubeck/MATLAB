function [uu,res] = sor(eqn, green, ps, F, b1, b2, h)
tol=1e-8;
its=100;

% Succesive Over-Relaxation
ub=ps(b1,b2);
uu=ub;
um=h*green(eqn(uu,F));
uu=uu+um;
normb=norm(eqn(ub,F),'fro');
i=1; res=norm(eqn(uu,F),'fro')/normb;
while i<its && res>=tol
    um=um+h*green(eqn(um,0));
    uu=uu+um;
    res=norm(eqn(uu,F),'fro')/normb;
    i=i+1;
end
end