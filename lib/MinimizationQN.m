classdef MinimizationQN < MinimizationND
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
    
    function self = MinimizationQN( fun, ls )
      self@MinimizationND( fun, ls ) ;
	    self.direction_short = 1e-3 ;
      self.angle_too_small = cos(pi/2+pi/180) ; % 90-1 degree
      self.method_names    = { ...
        'BFGS' , 'PSB' , 'DFP'
      } ;
      self.method = self.method_names{1} ;
    end
    
    function H = BFGS(self,g0,g1,d,alpha,H)
      y = g1 - g0;
      z = (H*y) / (d'*y);
      b = (alpha + y'*z) / (d'*y);
      H = H - (z*d' + d*z') + b*(d*d');
    end
    
    function H = PSB(self,g0,g1,d,alpha,H) % Powell Symmetric Broyden
      w = g1 + (alpha - 1)*g0;
      B = inv(H)+(d*w'+w*d')/(alpha*d'*d)-(d'*w)/(alpha)*(d*d');
      H = inv(B);
    end
    
    function H = DFP(self,g0,g1,d,alpha,H) % Davidon Fletcher and Powell rank 2 update
      y = g1 - g0;
      z = H*y;
      H = H + alpha*(d*d')/(d'*y) - (z*z')/(y'*z);  
    end
    
    
    function H = H_update(self,g0,g1,d,alpha,H)
      switch ( self.method )
        case 'BFGS' ; H = self.BFGS( g0,g1,d,alpha,H ) ;
        case 'PSB'  ;  H = self.PSB( g0,g1,d,alpha,H ) ;
        case 'DFP'  ; H = self.DFP( g0,g1,d,alpha,H ) ;
      end
    end
    
    function n = numOfMethods( self )
      % return the number of methods available
      n = length(self.method_names);
    end

    function name = activeMethod( self )
      % return active CG method used in minimization
      name = self.method;
    end

    function selectByNumber( self, k )
      % select the CG method by number
      if k < 1 || k > length(self.method_names)
        error('MinimizationCG, selectByNumber, k=%d out of range',k) ;
      end
      self.method = self.method_names{k} ;
    end

    function selectByName( self, name )
      % select the CG method by its name
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
      H         = inv(self.funND.hessian( xs ).') ;
      g1         = self.funND.grad( xs ).' ;      
      alpha     = 1 ;
      converged = false ;
      
      if self.debug_state
        self.x_history = reshape( x0, length(x0), 1 ) ;
      end
      
      for iter=1:self.max_iter
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
        % find search direction
        d = - H*g1;
        %
		% check direction ---------------------------------------------------------------
		if iter~=1 && abs(dot(g0,g1)) >= 0.2 * dot(g1,g1) % check restart criteria of Powell
          if self.debug_state
            fprintf(1,'Powell restart criteria, reset direction search, ...' ) ;
          end
          d = -g1 ; % reset direction as if H = eye
        else
          % check if the direction is descending, if direction is >= 89 degree from gradient reset to gradient directions
          % use >= to catch d == 0
          nrm_d  = norm(d) ;
          nrm_g1 = norm(g1) ;
          if dot( d, -g1 ) <= nrm_d * nrm_g1 * self.angle_too_small
            if self.debug_state
              fprintf(1,'direction angle about 90 degree, reset direction search, ...' ) ;
            end
            d = -g1 ;  % reset direction as if H = eye
          elseif nrm_d <= self.direction_short * nrm_g1
            if self.debug_state
              fprintf(1,'direction length too short, reset direction search, ...' ) ;
            end
            d = -g1 ; % reset direction as if H = eye
          end
        end
        % last check of direction search
        if norm(d,inf) == 0
          error('MinimizationCG, bad direction d == 0\n') ;
        end
		%--------------------------------------------------------------------------------------
        % minimize along search direction
        [xs,alpha] = self.step1D( xs, d, alpha ) ;
        %
        % update H (slied lesson 6 pag 61)
        g0 = g1;
        g1 = self.funND.grad( xs ).' ;
        H = self.H_update(g0,g1,d,alpha,H);
        
%         y = g1 - g0;
%         z = (H*y) / (d'*y);
%         b = (alpha + y'*z) / (d'*y);
%         H = H - (z*d' + d*z') + b*(d*d');
        
        if self.debug_state
          fprintf(1,' alpha = %8.4g\n', alpha) ;
        end
      end
    end

  end
end
