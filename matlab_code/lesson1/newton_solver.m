%
% declare the new function "newton_solver"
% that implement Newton method.
%
% INPUT Arguments:
%
% fun      = "pointer" to the function f(x)
% Dfun     = "pointer" to the function f`(x) (the derivative)
% x0       = starting approximation of alpha [f(alpha) = 0]
% tol      = tolerance to stopping criteria |f(x)| <= TOL -> STOP
% max_iter = maximum number of iterations
% verbose  = 'none' --> no message
%            'iter' --> print a message at each iteration
%
% INPUT Values:
%
%  x    = the computed approximation
%  iter = number of iteration done
%  flag = true  -> converged
%         false -> do not converge
%
function [x,iter,flag] = newton_solver( fun, Dfun, x0, tol, max_iter, verbose )
  x    = x0;    % initial approximation
  flag = false; % by default no convergence
  for iter=1:max_iter
    fx  = feval( fun, x ); % compute f(x);
    Dfx = feval( Dfun, x ); % compute f`(x);
    if Dfx == 0 % avoid division by 0
      if verbose == 'iter'
        fprintf( 'find Df(x) == 0\n');
      end
      flag = false;
      return;
    end
    % perform Newton step
    x = x - fx/Dfx;
    % print iteration
    if verbose == 'iter'
      fprintf( 'iter=%2d x=%g f(x)=%g\n', iter, x, fx );
    end
    % check convergence
    if abs(fx) <= tol
      if verbose == 'iter'
        fprintf('convergence reached at iter %d\n',iter);
      end
      flag = true;
      return;
    end
  end
end
