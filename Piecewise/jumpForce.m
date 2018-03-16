function [f] = jumpForce(xi,xe,x0,x1,A0,jumps)
% Computes the right hand side for a dirac Delta function 
if any(xe==xi)
    f=(x1==xi);
else
    p=size(A0,1);
    d1=1:p;
    d2=p+d1;
    % Locate singularity
    [~,imax]=max(xe>xi);
    j=imax-1;
    jd=(1+(j-1)*(p-1)):(1+j*(p-1));
    
    % Subdomain
    xa=x1(jd(1));
    xb=x1(jd(end));
    % Jacobians
    h0=(xb-xa)/2;
    h1=(xi-xa)/2;
    h2=(xb-xi)/2;
    % Centers
    c1=(xi+xa)/2;
    c2=(xi+xb)/2;
    % Supergrid
    xs=[c1+h1*x0; c2+h2*x0];
    % Canonical reference frame [-1,1]
    xx=2/(xb-xa)*(xs-(xb+xa)/2);
    xi0=2/(xb-xa)*(xi-(xb+xa)/2);
    % Compute piecewise corrections
    jumps=jumps.*(h0.^(-1:numel(jumps)-2))';
    [s1,s2]=piecewiseLagrange(x0,xi0,jumps);
    % Interpolator operator
    C=legC(x0,xx);
    % Compute force
    fdisc=-C'*[(h0/h1)*A0*C(d1,:)*s1; (h0/h2)*A0*C(d2,:)*s2];
    f0=fdisc+jumps(2)*C(p+1,:)';
    f=zeros(size(x1));
    f(jd)=f0;
end
end

