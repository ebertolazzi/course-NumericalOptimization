%
% solve the problem f(x) = 0 starting from point x0
%
function x1 = solve( self, f, x0, x1 )
  self.setup( f );
  self.flag = false;
  self.x_history = [ x0 x1 ]; % save iteration
  fx0 = feval( self.fun, x0 ); % compute f(x);
  fx1 = feval( self.fun, x1 ); % compute f(x);
  for iter=1:self.max_iter
    % check convergence
    if abs(fx1) <= self.tol
      if self.verbose == 'iter'
        fprintf('convergence reached at iter %d\n',iter);
      end
      self.flag = true;
      break;
    end
    top = fx1;
    bot = (fx1-fx0)/(x1-x0);
    if bot == 0 % avoid division by 0
      if self.verbose == 'iter'
        fprintf( 'find f`(x) == 0\n');
      end
      break;
    end
    % perform Secant step
    x2  = x1 - top/bot;
    fx2 = feval( self.fun, x2 ); % compute f(x);
    % print iteration
    if self.verbose == 'iter'
      fprintf( 'iter=%2d x=%g f(x)=%g\n', iter, x2, fx2 );
    end
    % save iteration
    self.x_history = [self.x_history x2];
    % shift values
    x0  = x1;   x1 = x2;
    fx0 = fx1; fx1 = fx2;
  end
  self.iter = iter;
end
