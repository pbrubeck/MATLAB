function y = count_sign_changes(x,tol)

%COUNT_SIGN_CHANGES counts the number of sign changes
%
% y = COUNT_SIGN_CHANGES(x,tol) returns the number of changes of sign
% in real part of vector x
% components with absolute value less then tol (1e-6) do not count

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 8.9.2015

if nargin<2, tol=1e-6; end

if norm(real(x))>norm(imag(x))
    z = real(x)/norm(real(x),'inf');
else
    z = imag(x)/norm(imag(x),'inf');
end
z = sign(z).*(abs(z)>tol);
z = z(z~=0);
y = sum(abs(z(2:end)-z(1:end-1)))/2;
    

 