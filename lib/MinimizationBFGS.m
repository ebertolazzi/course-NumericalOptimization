classdef MinimizationBFGS < MinimizationND
  %
  % ADD comment describing the class what do etc )
  %
  %
  %

  methods
    function self = MinimizationBFGS( fun, ls )
      self@MinimizationND( fun, ls ) ;
    end
    
    function [xs,converged] = minimize( self, x0 )
      xs        = x0 ;
      H         = inv(self.funND.hessian( xs ).') ;
      H         = eye(length(x0));
      g         = self.funND.grad( xs ).' ;      
      alpha     = 1 ;
      converged = false ;
      
      if self.debug_state
        self.x_history = reshape( x0, length(x0), 1 ) ;
      end
      
      for iter=1:self.max_iter

        
        % find search direction
        d = - H*g;
        
        % check if converged
        nrm = norm(d,inf) ;
        converged = nrm < self.tol ;
        if converged
          if self.debug_state
            fprintf(1,'solution found, ||grad f||_inf = %g < %g\n', nrm, self.tol ) ;
          end
          break ;
        end 
        %
        % only for debug
        if self.debug_state
          fprintf(1,'iter = %5d ||grad f||_inf = %12.6g ...', iter, nrm ) ;
        end
        %
        % minimize along search direction
        [xs,alpha] = self.step1D( xs, d, alpha ) ;
        %
        % update H (slied lesson 6 pag 61)
        grad_up = self.funND.grad( xs ).' ;
        y = grad_up - g;
        z = (H*y) / (d'*y);
        g = grad_up;
        b = (alpha + y'*z) / (d'*y);
        H = H - (z*d' + d*z') + b*(d*d');
        
        if self.debug_state
          fprintf(1,' alpha = %8.4g\n', alpha) ;
        end
      end
    end

  end
end
