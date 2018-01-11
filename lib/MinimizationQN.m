classdef MinimizationQN < MinimizationND
 %
  % Description
  % -----------
  % Minimization of a nonlinear multiuvariate functions using QUASI NEWTON methods.
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
    method          % selected QN method
    method_names    % list of available methods
    direction_short %
    angle_too_small %
  end
  
   
  methods (Hidden = true)
    
	
	function d = BFGS( self, g1, g0, H0, d0, s, alpha )
      % References
      % (slides lesson 6 pag 61)
      y = g1 - g0;
      z = (H0*y) / (d0'*y);
      b = (alpha + y'*z) / (d0'*y);
      H1 = H0 - (z*d0' + d0*z') + b*(d0*d0');
	  d  = H1*g1;
    end
	
	
	function d = dirEval( self, g1, g0, H0, d0, s, alpha )
      % Compute the value of `beta` for the selected QUASI NEWTON method
      switch ( self.method )
      case 'BFGS'  ; d = self.BFGS( g1, g0, H0, d0, s, alpha )
      otherwise
        error('MinimizationCG, dirEval no method selected (this should not be happen)');
      end
    end
	
  end
	
	
  methods	
    
    function self = MinimizationQN( fun, ls )
	  %
      % fun = function to be minimized
      % ls  = linesearch used, can be LinesearchArmijo, LinesearchWolfe, LinesearchGoldenSection
      self@MinimizationND( fun, ls ) ;
	  self.method_names    = { ...
        'BFGS'...
	  };
	  self.method          = self.method_names{1} ;
	  self.direction_short = 1e-3 ;
      self.angle_too_small = cos(pi/2+pi/180) ; % 90-1 degree
    end

	function n = numOfMethods( self )
      % return the number of methods available
      n = length(self.method_names);
    end 
	  
	  
	function name = activeMethod( self )
      % return active QN method used in minimization
      name = self.method;
    end 

	
	function selectByNumber( self, k )
      % select the QN method by number
      if k < 1 || k > length(self.method_names)
        error('MinimizationCG, selectByNumber, k=%d out of range',k) ;
      end
      self.method = self.method_names{k} ;
    end	
	
    function selectByName( self, name )
      % select the QN method by its name
      if ischar(name)
        for k=1:length(self.method_names)
          if strcmp(name,self.method_names{k})
            self.method = name ;
            return ;
          end
        end
        error('MinimizationCG, selectByName, name=%s not found',name) ;
      else
        error('MinimizationCG, selectByName, expected string as arument, found %s',class(name)) ;
      end
    end	
	 
    function [xs,converged] = minimize( self, x0 )
	  xs        = x0 ;
	  alpha     = 1 ;
      converged = false ;
      if self.debug_state
        self.x_history = reshape( x0, length(x0), 1 ) ;
      end
      
	  
      for iter=1:self.max_iter
        %
        % gradient of the function
        g1 = self.funND.grad( xs ).' ;
		
		% inverse of the hessian of the function
		H1 = inv(self.funND.hessian( xs ).') ;
        %
        % check if converged
        nrm_g1    = norm(g1,inf) ;
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
          H1 = eye(self.funND.N);
		  d = - H1*g1 ; % first iteration, search direction is -hessian^-1*gradient
        elseif abs(dot(g0,g1)) >= 0.2 * dot(g1,g1) % check restart criteria of Powell
          if self.debug_state
            fprintf(1,'Powell restart criteria, reset direction search, ...' ) ;
          end
          d = - H1*g1 ; % reset direction    
        else
		  d = self.dirEval( g1, g0, H0, d, s, alpha ) ; % s = x1 - x0 last step! -----------------------------------------------------
          % check if the direction is descending, if direction is >= 89 degree from hessian^(-1)*gradient reset to gradient previous
          % use >= to catch d == 0
          nrm_d  = norm(d) ;
          nrm_d_old = norm(H1*g1) ;
          if dot( d, -H1*g1 ) <= nrm_d * nrm_d_old * self.angle_too_small
            if self.debug_state
              fprintf(1,'direction angle about 90 degree, reset direction search, ...' ) ;
            end
            d = - H1*g1 ; % reset direction
          elseif nrm_d <= self.direction_short * nrm_d_old
            if self.debug_state
              fprintf(1,'direction length too short, reset direction search, ...' ) ;
            end
            d = - H1*g1 ; % reset direction
          end
        end
		%
        % minimize along search direction
        dot(H1*g1,g1)
        [xs,alpha,ok] = self.step1D( xs, d, alpha ) ;
        if ~ok
          % step failed try to use gradient direction ---------------------------------
          d = - g1;
          [xs,alpha,ok] = self.step1D( xs, d, alpha ) ;
          if ~ok
            % cannot advance see if accept a low precision solution
            warning('MinimizationCG, step1D failed') ;
            return ;
          end
        end
        %
	    if self.debug_state
          fprintf(1,' alpha = %8.4g\n', alpha) ;
        end	
        %
        % save old gradient  and hessian^(-1) and step
        g0 = g1 ;
		H0 = H1 ;
        s  = alpha*d ;
		  
      end
    end

  end
end
