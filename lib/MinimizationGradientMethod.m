classdef MinimizationGradientMethod < MinimizationND

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = MinimizationGradientMethod( fun, ls )
      self@MinimizationND( fun, ls );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, converged ] = minimize( self, x0 )
      x     = x0(:);
      alpha = 1;
      converged = false;
      if self.debug_state
        self.x_history = x(:);
      end
      for iter=1:self.max_iter
        %
        % find search direction = - gradient of the function
        d = -self.funND.grad( x );
        %
        % check if converged
        norm_inf_d = norm(d,inf);
        converged  = norm_inf_d < self.tol;
        if converged
          if self.debug_state
            fprintf('solution found, ||grad f||_inf = %g < %g\n', ...
                     norm_inf_d, self.tol );
          end
          break;
        end
        %
        % only for debug
        if self.debug_state
          fprintf('iter = %5d ||grad f||_inf = %12.6g ...', ...
                   iter, norm_inf_d );
        end
        %
        % minimize along search direction
        [ x, alpha ] = self.step1D( x, d(:), 10*alpha );
        %
        if self.debug_state
          fprintf(' alpha = %8.4g\n', alpha);
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
