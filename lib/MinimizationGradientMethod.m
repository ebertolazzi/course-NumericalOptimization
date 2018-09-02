classdef MinimizationGradientMethod < MinimizationND

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = MinimizationGradientMethod( fun, ls )
      self@MinimizationND( fun, ls );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [xs,converged] = minimize( self, x0 )
      xs    = x0;
      alpha = 1;
      converged = false;
      if self.debug_state
        self.x_history = reshape( x0, length(x0), 1 );
      end
      for iter=1:self.max_iter
        %
        % find search direction = - gradient of the function
        d = -self.funND.grad( xs ).';
        %
        % check if converged
        norm_inf_d = norm(d,inf);
        converged  = norm_inf_d < self.tol;
        if converged
          if self.debug_state
            fprintf(1,'solution found, ||grad f||_inf = %g < %g\n', norm_inf_d, self.tol );
          end
          break;
        end
        %
        % only for debug
        if self.debug_state
          fprintf(1,'iter = %5d ||grad f||_inf = %12.6g ...', iter, norm_inf_d );
        end
        %
        % minimize along search direction
        [xs,alpha] = self.step1D( xs, d, alpha );
        %
        if self.debug_state
          fprintf(1,' alpha = %8.4g\n', alpha);
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
