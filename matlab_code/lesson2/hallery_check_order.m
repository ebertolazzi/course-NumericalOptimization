fun      = @(x) x.*cos(x) + x + x./exp(x) - cos(x) - 1 - 1./exp(x);
Dfun     = @(x) cos(x) - x.*sin(x) + 1 + 2./exp(x) - x./exp(x) + sin(x);
DDfun    = @(x)-2*sin(x) - x*cos(x) - 3./exp(x) + x./exp(x) + cos(x);
x0       = -3;
HALLEY = halley_class();
HALLEY.set_tolerance( 1e-10 );
x = HALLEY.solve( fun, Dfun, DDfun, x0 );
fprintf('x = %g\n',x);
fprintf('f(x) = %g\n',fun(x));

xh = HALLEY.get_history();
% compute error
sol = 1;
err = abs(xh-sol);
% log(e(k+2)/e(k+1))/log(e(k+1)/e(k))
p = log(err(3:end)./err(2:end-1))./log(err(2:end-1)./err(1:end-2));
p.'
