classdef MinimizationQuasiNewton < MinimizationND
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
    method
    method_names
    direction_short %
    angle_too_small %
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = MinimizationQuasiNewton( fun, ls )
      self@MinimizationND( fun, ls );
	  self.direction_short = 1e-10;
      self.angle_too_small = cos(pi/2+pi/180); % 90-1 degree
      self.method_names = { 'BFGS', 'DFP' };
      self.method       = self.method_names{1};
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function H = H_update( self, g0, g1, s, H )
      y  = g1 - g0;
      sy = dot(s,y);
      if sy > 0
        z = H*y;
        switch ( self.method )
        case 'BFGS';
          H = H + (dot(y,s+z)/sy^2)*(s*s') - (z*s' + s*z')/sy;
        case 'DFP';
          H = H - (z/dot(z,y))*z'+ (s/sy)*s';
        otherwise;
          error( 'MinimizationQN, method `%s` not supported', self.method );
        end
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function n = numOfMethods( self )
      % return the number of methods available
      n = length(self.method_names);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function name = activeMethod( self )
      % return active CG method used in minimization
      name = self.method;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function selectByNumber( self, k )
      % select the CG method by number
      if k < 1 || k > length(self.method_names)
        error('MinimizationCG, selectByNumber, k=%d out of range',k);
      end
      self.method = self.method_names{k};
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function selectByName( self, name )
      % select the CG method by its name
      if ischar(name)
        for k=1:length(self.method_names)
          if strcmp(name,self.method_names{k})
            self.method = name;
            return;
          end
        end
        error('MinimizationCG, selectByName, name=%s not found',name);
      else
        error('MinimizationCG, selectByName, expected string as arument, found %s',class(name));
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x, converged ] = minimize( self, x0 )
      self.n_fail = 0;
      x           = x0(:);
      H           = eye(length(x));
      alpha       = 1;
      converged   = false;

      if self.save_iterate
        self.x_history = x;
      end

      for iter=1:self.max_iter
        self.iter = iter;
        g1 = self.funND.grad( x ).';
        % check if converged
        norm_inf_g1 = norm(g1,inf);
        converged   = norm_inf_g1 < self.tol;
        if converged
          if self.verbose
            fprintf( 'solution found, ||grad f||_inf = %g < %g\n', ...
                     norm_inf_g1, self.tol );
          end
          break;
        end
        %
        if self.verbose
          fprintf( '[%s] iter = %5d ||grad f||_inf = %12.6g', ...
                   self.method, self.iter, norm_inf_g1 );
        end
        %
        if self.iter > 1
          H = self.H_update( g0, g1, s, H );
        end
        %
        % find search direction
        d = -H*g1;
        %
        % check direction -------------------------------------------------
        if self.iter > 1 && abs(dot(g0,g1)) >= 0.2 * dot(g1,g1) % check restart criteria of Powell
          fprintf(2,' [Powell restart] ');
          d = -g1; % reset direction as if H = eye
        else
          % check if the direction is descending,
          % if direction is >= 89 degree from gradient reset to gradient directions
          % use >= to catch d == 0
          norm_2_d  = norm(d);
          norm_2_g1 = norm(g1);
          if dot( d, -g1 ) <= norm_2_d * norm_2_g1 * self.angle_too_small
            fprintf(2,' [angle ~90, reset] ');
            d = -g1; % reset direction as if H = eye
          end
        end
        % last check of direction search
        if norm(d,inf) == 0
          error('MinimizationQN, bad direction d == 0\n');
        end
        %------------------------------------------------------------------
        % minimize along search direction
        [ x1, alpha, ok ] = self.step1D( x, d, 10*alpha );
        if ~ok
          % cannot advance see if accept a low precision solution
          fprintf(2,'\nMinimizationQN, step1D failed\n');
          return;
        end
        %
        fprintf( ' alpha = %8.4g\n', alpha );
        s  = x1-x;
        x  = x1;
        g0 = g1;
      end
    end
  end
end
