%
% solve the problem f(x) = 0 starting from point x0
%
function x = solve( self, f, Df, DDf, x0 )
  self.setup(f,Df,DDf);
  self.flag = false;
  x = x0;    % initial approximation
  for iter=1:self.max_iter
    fx   = feval( self.fun,   x ); % compute f(x);
    Dfx  = feval( self.Dfun,  x ); % compute f`(x);
    DDfx = feval( self.DDfun, x ); % compute f``(x);
    top = 2*fx*Dfx;
    bot = 2*Dfx^2-fx*DDfx;
    if bot == 0 % avoid division by 0
      if self.verbose == 'iter'
        fprintf( 'find 2*f`(x)^2-f(x)*f``(x) == 0\n');
      end
      break;
    end
    % perform Halley step
    x = x - top/bot;
    % print iteration
    if self.verbose == 'iter'
      fprintf( 'iter=%2d x=%g f(x)=%g\n', iter, x, fx );
    end
    % check convergence
    if abs(fx) <= self.tol
      if self.verbose == 'iter'
        fprintf('convergence reached at iter %d\n',iter);
      end
      self.flag = true;
      break;
    end
  end
  self.iter   = iter;
  self.x_last = x;
end
