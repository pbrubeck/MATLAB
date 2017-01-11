function L = LapBel(f,m)

% LAPBEL Laplace-Beltrami operator.
%    L = LAPBEL(f,m) creates an anonymous function L for evaluating the
%    Laplace-Beltrami operator of the function f defined on the manifold m,
%    where both f and m are given as anonymous functions.
%
%    The manifold m can be specified either in parametrized or implicit form.
%    In the first case, to define an n-dimensional manifold in R^d, m is a function
%    mapping the parameters u1,...,un to a cell array of d coordinate functions.
%    In the second case, m is a real-valued function depending on u1,...,un
%    defining a hypersurface in R^n by the condition m(u1,...,un) = const.
%    In both cases, f is a function depending on u1,...,un.
%
%    LAPBEL is based on the <a href="http://www.mathworks.com/matlabcentral/fileexchange/56856-autodiff">AutoDiff</a> toolbox.
%
%    Example: Torus colored by Laplace-Beltrami of some function.
%    m = @(u,v) {(3+cos(u)).*cos(v), (3+cos(u)).*sin(v), sin(u)};
%    f = @(u,v) sin(3*v-2*u);
%    L = LapBel(f,m);
%    [u,v] = ndgrid(linspace(-pi,pi,200));
%    M = m(u,v);
%    figure, surf(M{1},M{2},M{3},L(u,v))
%    axis equal, grid on, shading interp, camlight, colormap(prism(55))
%
%    Example: Ellipsoid colored by Laplace-Beltrami of some function.
%    m = @(x,y,z) 2*x.^2 + y.^2 + 4*z.^2;
%    f = @(x,y,z) sqrt(x.^2 + y.^2 + z.^2);
%    [x,y,z] = meshgrid(linspace(0,.71),linspace(0,1),linspace(0,.5));
%    L = LapBel(f,m);
%    figure, isosurface(x,y,z,m(x,y,z),1,L(x,y,z),'noshare')
%    colormap(jet(11)), axis equal, grid on, camlight, view(130,30)

%    Ulrich Reif
%    May 31, 2016

if ~isa(f,'function_handle')
  if f.n==2                      % parametric
    m = [m{:}];
    J = ajac(m(:));
    G = J'*J;
    g = sqrt(det(G));
    L = aeval(adiv(G\agrad(f)*g)/g,0);
  else                           % implicit
    N = agrad(m);
    N = N/norm(N);
    L = aeval(alap(f) - N'*(ahess(f)*N + agrad(f)*adiv(N)),0);
  end
else
  fs = func2str(f);
  d = true(1,length(find(fs(1:find(fs==')',1))==','))+1);
  [a1,a2,p] = prepfct(ainit([]),d,2);
  eval(['L = @(' a1 ') LapBel(f(' a2 '),m(' a2 '));'])
end
