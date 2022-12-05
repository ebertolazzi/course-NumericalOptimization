fun      = @(x) x.*cos(x) + x + x./exp(x) - cos(x) - 1 - 1./exp(x);
x0       = -2;
x1       = -1;
%tol      = 1e-50;
%max_iter = 100;
%verbose  = 'iter';
SE = secant_class();
SE.set_tolerance( 1e-10 );
x = SE.solve( fun, x0, x1 );
fprintf('x = %g\n',x);
fprintf('f(x) = %g\n',fun(x));

xh = SE.get_history();
% compute error
sol = 1;
err = abs(xh-sol);
% log(e(k+2)/e(k+1))/log(e(k+1)/e(k))
p = log(err(3:end)./err(2:end-1))./log(err(2:end-1)./err(1:end-2));
p.'
