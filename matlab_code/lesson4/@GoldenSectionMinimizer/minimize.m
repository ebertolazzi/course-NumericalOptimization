function [a,b] = minimize( self, fun, a_in, b_in )
  % reorder a and b in such a way b > a
  a = min(a_in,b_in);
  b = max(a_in,b_in);
  self.history = [a;b];
  % early exit, check tolerance
  if b-a <= self.tol
    return; % the length of [a,b] is less than tolerance, exiting from the loop
  end
  %
  %   build the sub intervals [ a mu ] and [ lambda b ]
  %
  %                c = (a+b)/2
  %   |--------+---|---+--------|
  %   a     lambda     mu       b
  %
  % compute the length of the subintervals [ golden section * length ]
  dlen = self.tau*(b-a);
  % compute egdes of sub intervals
  lambda = b-dlen;
  mu     = a+dlen;
  % compute values at `lambda` and `mu`
  fl     = feval(fun,lambda);
  fm     = feval(fun,mu);
  for k=1:self.max_iter
    % new length of the next interval after selection
    dlen = dlen*self.tau;
    if fl > fm
      % Select interval [lambda,b] --> fl > fm <= fb
      %
      % the new interval is [a,b] <- [lambda,b]
      a = lambda;
      %
      % Compute the next lambda and mu with f(lambda) and f(mu)
      % new_lambda = old_mu
      %
      lambda = mu;
      fl     = fm;
      %
      % new_mu  must be computed
      %
      mu = a+dlen;
      fm = feval(fun,mu);
    else
      % Select interval [a,mu] --> fa >= fl <= fm
      %
      % the new interval is [a,b] <- [a,mu]
      b = mu;
      %
      % Compute the next lambda and mu with f(lambda) and f(mu)
      % new_mu = old_lambda
      %
      mu = lambda;
      fm = fl;
      %
      % new_lambda must be computed
      %
      lambda = b-dlen;
      fl     = feval(fun,lambda);
    end
    self.history = [self.history [a;b]];
    if b-a <= self.tol
      break; % the length of [a,b] is less than tolerance, exiting from the loop
    end
  end
end
