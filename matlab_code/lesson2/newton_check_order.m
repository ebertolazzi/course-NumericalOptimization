fun      = @(x) x.*cos(x) + x + x./exp(x) - cos(x) - 1 - 1./exp(x);
Dfun     = @(x) cos(x) - x.*sin(x) + 1 + 2./exp(x) - x./exp(x) + sin(x);
x0       = -1;
%tol      = 1e-50;
%max_iter = 100;
%verbose  = 'iter';
NS = newton_class();
NS.set_tolerance( 1e-10 );
x = NS.solve( fun, Dfun, x0 );
fprintf('x = %g\n',x);
fprintf('f(x) = %g\n',fun(x));

xh = NS.get_history();
% compute error
sol = 1;
err = abs(xh-sol);
% log(e(k+2)/e(k+1))/log(e(k+1)/e(k))
p = log(err(3:end)./err(2:end-1))./log(err(2:end-1)./err(1:end-2));
p.'