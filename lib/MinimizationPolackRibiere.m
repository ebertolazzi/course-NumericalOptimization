classdef MinimizationPolackRibiere < MinimizationND

  methods
    function self = MinimizationPolackRibiere( fun, ls )
      self@MinimizationND( fun, ls ) ;
    end

    %
    % beta depends on actual and previous gradient
    % beta depends also on previous direction
    %
    function beta = betaEval( self, g1, g0, d0 )
      % Polack Ribiere (do not use d0)
      beta = max(0,dot(g1,g1-g0)/dot(g0,g0)) ;
      % Hestenes Stiefel
      %beta = dot(g1,g1-g0)/dot(g1-g0,d0) ;
    end

    function [xs,converged] = minimize( self, x0 )
      xs    = x0 ;
      alpha = 1 ;
      converged = false ;
      if self.debug_state
        self.x_history = reshape( x0, length(x0), 1 ) ;
      end
      %
      for iter=1:self.max_iter
        %
        % gradient of the function
        g1 = self.funND.grad( xs ).' ;
        %
        % check if converged
        nrm_g1 = norm(g1,inf) ;
        converged = nrm_g1 < self.tol ;
        if converged
          if self.debug_state
            fprintf(1,'solution found, ||grad f||_inf = %g < %g\n', nrm_g1, self.tol ) ;
          end
          break ;
        end 
        %
        % only for debug
        if self.debug_state
          fprintf(1,'iter = %5d ||grad f||_inf = %12.6g ...', iter, nrm_g1 ) ;
        end
        %
        % build search direction
        %
        if iter == 1
          d1 = -g1 ; % first iteration, search direction is -gradient
        else
          beta = self.betaEval( g1, g0, d0 ) ;
          d1 = -g1 + beta*d0 ;
          % check if the direction is descending
          if dot( d1, g1 ) > -norm(d1) * nrm_g1 * sqrt(eps)
            if self.debug_state
              fprintf(1,'reset direction search, ...' ) ;
            end
            d1 = -g1 ; % reset direction
          end
        end
        %
        % minimize along search direction
        [xs,alpha] = self.step1D( xs, d1, alpha ) ;
        %
        if self.debug_state
          fprintf(1,' alpha = %8.4g\n', alpha) ;
        end
        %
        % save old gradient and direction
        g0 = g1 ;
        d0 = d1 ;
      end
    end

  end
end
