%
% solve the problem f(x) = 0 starting from point x0
%
function x = solve( self, f, Df, x0 )
  self.setup( f, Df );
  self.flag = false;
  x = x0;              % initial approximation
  self.x_history = x0; % save iteration
  for iter=1:self.max_iter
    fx = feval( self.fun, x ); % compute f(x);
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
    Dfx = feval( self.Dfun,  x ); % compute f`(x);
    top = fx;
    bot = Dfx;
    if bot == 0 % avoid division by 0
      if self.verbose == 'iter'
        fprintf( 'find f`(x) == 0\n');
      end
      break;
    end
    % perform Newton step
    x = x - top/bot;
    % save iteration
    self.x_history = [self.x_history x];
  end
  self.iter = iter;
end
