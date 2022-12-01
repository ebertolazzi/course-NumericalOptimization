fun      = @(x) x.*exp(-x.*x)-sin(x);
Dfun     = @(x) exp(-x.^2) - 2.*x.^2*exp(-x.^2) - cos(x);
x0       = 1;
tol      = 1e-50;
max_iter = 100;
verbose  = 'iter';
[x,iter,flag] = newton_solver( fun, Dfun, x0, tol, max_iter, verbose );
fprintf('x = %g\n',x);