classdef MinimizationBFGS < MinimizationND
 %
  % Description
  % -----------
  % Minimization of a nonlinear multiuvariate functions using (BFGS) method.
  % The algorithm is describbed in the references.
  %
  % References
  % ----------
  %
  % @article{2006,
  %     author   = {},
  %     journal  = {},
  %     keywords = {},
  %     number   = {},
  %     pages    = {},
  %     title    = {},
  %     volume   = {},
  %     year     = {}
  % }
  % @article{2006,
  %     author  = {},
  %     title   = {},
  %     journal = {},
  %     volume  = {},
  %     number  = {},
  %     year    = {},
  %     pages   = {},
  %     doi     = {},
  % }
  %
  % Authors: Enrico Bertolazzi & Davide Vignotto
  %
  properties (SetAccess = private, Hidden = true)
    angle_too_small
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = MinimizationBFGS( fun, ls )
      self@MinimizationND( fun, ls );
      self.angle_too_small = cos(pi/2+pi/180); % 90-1 degree
    end

    function [ x,converged] = minimize( self, x0 )
      x         = x0(:);
      %H         = inv(self.funND.hessian(x));
      H         = eye(length(x));
      alpha     = 1;
      converged = false;

      if self.debug_state
        self.x_history = x;
      end

      for iter=1:self.max_iter
        g1 = self.funND.grad( x ).';
        % check if converged
        norm_inf_g1 = norm(g1,inf);
        converged   = norm_inf_g1 < self.tol;
        if converged
          if self.debug_state
            fprintf(1,'solution found, ||grad f||_inf = %g < %g\n', norm_inf_g1, self.tol );
          end
          break;
        end
        %
        % only for debug
        if self.debug_state
          fprintf(1,'iter = %5d ||grad f||_inf = %12.6g ...', iter, norm_inf_g1 );
        end
        %
        if iter > 1
          y  = g1 - g0;
          sy = dot(s,y);
          if sy > 0
            z = H*y;
            if true
              % BFGS update
              H = H + (dot(y,s+z)/sy^2)*s*s' - (z*s' + s*z')/sy;
            else
              % DFP update
              H = H - (z/dot(z,y))*z'+ (s/sy)*s';
            end
          end
        end
		%
        % find search direction
        d = -H*g1;
        %
		% check direction ---------------------------------------------------------------
		if iter~=1 && abs(dot(g0,g1)) >= 0.2 * dot(g1,g1) % check restart criteria of Powell
          if self.debug_state
            fprintf(1,'Powell restart criteria, reset direction search, ...' );
          end
          d = -g1; % reset direction as if H = eye
        else
          % check if the direction is descending,
          % if direction is >= 89 degree from gradient reset to gradient directions
          % use >= to catch d == 0
          norm_2_d  = norm(d);
          norm_2_g1 = norm(g1);
          if dot( d, -g1 ) <= norm_2_d * norm_2_g1 * self.angle_too_small
            if self.debug_state
              fprintf( ' [ angle ~90 degree, reset] ' );
            end
            d = -g1; % reset direction as if H = eye
          end
        end
        % last check of direction search
        if norm(d,inf) == 0
          error('MinimizationCG, bad direction d == 0\n');
        end
		%--------------------------------------------------------------------------------------
        % minimize along search direction
        [ x, alpha, ok ] = self.step1D( x, d, 10*alpha );
        if ~ok
          % step failed try to use gradient direction
          d = -g1;
          H = eye(length(d));
          [x,alpha,ok] = self.step1D( x, d, alpha );
          if ~ok
            % cannot advance see if accept a low precision solution
            warning('MinimizationBFGS, step1D failed');
            return;
          end
        end
        %
        if self.debug_state
          fprintf(1,' alpha = %8.4g\n', alpha);
        end
        s  = alpha*d;
        g0 = g1;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
