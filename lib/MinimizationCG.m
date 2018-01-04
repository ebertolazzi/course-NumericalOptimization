classdef MinimizationCG < MinimizationND
  %
  % Description
  % -----------
  % Minimization of a nonlinear multiuvariate functions using nonlinbear conjugate gradient mnethods.
  % The algorithm is describbed in the references.
  %
  % References
  % ----------
  %
  % @article{2006,
  %     author   = {Hager, William W. and Zhang, Hongchao},
  %     journal  = {Pacific journal of Optimization},
  %     keywords = {bibliot-11-08-07},
  %     number   = {1},
  %     pages    = {35--58},
  %     title    = {A survey of nonlinear conjugate gradient methods},
  %     volume   = {2},
  %     year     = {2006}
  % }
  % @article{2006,
  %     author  = {Hager, William W. and Zhang, Hongchao},
  %     title   = {Algorithm 851: CG_DESCENT, a Conjugate Gradient Method with Guaranteed Descent},
  %     journal = {ACM Trans. Math. Softw.},
  %     volume  = {32},
  %     number  = {1},
  %     year    = {2006},
  %     pages   = {113--137},
  %     doi     = {10.1145/1132973.1132979},
  % }
  %
  % Author: Enrico Bertolazzi
  %

  properties (SetAccess = private, Hidden = true)
    method          % selected nonliner CG method
    method_names    % list of available methods
    direction_short %
    angle_too_small %
  end

  methods (Hidden = true)

    function beta = HS( self, g1, g0, d0 )
      % M.R.Hestenes and E.L.Stiefel  
      % Methods of conjugate gradients for solving linear systems,
      % J. Research Nat. Bur. Standards, 49 (1952), pp. 409–436
      y    = g1 - g0 ;
      beta = dot(g1,y)/dot(d0,y) ;
    end

    function beta = HSplus( self, g1, g0, d0 )
      % M. J. D. Powell
      % Nonconvex minimization calculations and the conjugate gradient method,
      % Numerical Analysis (Dundee, 1983), Lecture Notes in Mathematics, Vol. 1066,
      % Springer-Verlag, Berlin, 1984, pp. 122–141.
      y    = g1 - g0 ;
      beta = max(0,dot(g1,y)/dot(d0,y)) ;
    end

    function beta = FR( self, g1, g0, d0 )
      % R. Fletcher and C. Reeves
      % Function minimization by conjugate gradients,
      % Comput. J., 7 (1964), pp. 149–154.
      beta = dot(g1,g1)/dot(g0,g0) ;
    end

    function beta = PRP( self, g1, g0, d0 )
      % E. Polak and G. Ribiere
      % Note sur la convergence de directions conjugees,
      % Rev. Francaise Informat Recherche Opertionelle, 3e Annee 16 (1969), pp. 35–43.
      %
      % B. T. Polyak,
      % The conjugate gradient method in extreme problems,
      % USSR Comp. Math. Math. Phys., 9 (1969), pp. 94–112.
      beta = dot(g1,g1-g0)/dot(g0,g0) ;
    end

    function beta = PRPplus( self, g1, g0, d0 )
      % J. C. Gilbert and J. Nocedal
      % Global convergence properties of conjugate gradient methods for optimization,
      % SIAM J. Optim., 2 (1992), pp. 21–42.
      beta = max(0,dot(g1,g1-g0)/dot(g0,g0)) ;
    end

    function beta = CD( self, g1, g0, d0 )
      % R. Fletcher
      % Practical Methods of Optimization vol. 1: Unconstrained Optimization,
      % John Wiley & Sons, New York, 1987.
      beta = -dot(g1,g1)/dot(d0,g0);
    end

    function beta = LS( self, g1, g0, d0 )
      % Y. Liu and C. Storey
      % Efficient generalized conjugate gradient algorithms, Part 1: Theory,
      % J. Optim. Theory Appl., 69 (1991), pp. 129–137.
      beta = -dot(g1,g1-g0)/dot(d0,g0);
    end

    function beta = DY( self, g1, g0, d0 )
      % Y. H. Dai and Y. Yuan
      % A nonlinear conjugate gradient method with a strong global convergence property
      % SIAM J. Optim., 10 (1999), pp. 177–182.
      beta = dot(g1,g1)/dot(d0,g1-g0) ;
    end

    function beta = N( self, g1, g0, d0 )
      % W. W. Hager and H. Zhang
      % A new conjugate gradient method with guaranteed descent and an efficient line search
      % SIAM J. Optim., 16(1), 170–192. 2003
      y    = g1 - g0 ;
      tmp  = dot(d0,y) ; 
      z    = y - (2*dot(y,y)/tmp)*d0 ;
      beta = dot(z,g1)/tmp ; 
    end

    function beta = HyTAS( self, g1, g0, d0 )
      % D. Touati-Ahmed and C. Storey
      % Efficient hybrid conjugate gradient techniques,
      % J. Optim. Theory Appl., 64 (1990), pp. 379–397.
      g0g0    = dot(g0,g0) ;
      betaPRP = dot(g1,g1-g0)/g0g0 ;
      betaFR  = dot(g1,g1)/g0g0 ;
      if 0 <= betaPRP && betaPRP <= betaFR
        beta = betaPRP ;
      else
        beta = betaFR ;
      end
    end

    function beta = HyHS( self, g1, g0, d0 )
      % Y. F. Hu and C. Storey
      % Global convergence result for conjugate gradient methods,
      % J. Optim. Theory Appl., 71 (1991), pp. 399–405.
      g0g0    = dot(g0,g0) ;
      betaPRP = dot(g1,g1-g0)/g0g0 ;
      betaFR  = dot(g1,g1)/g0g0 ;
      beta    = max(0,min(betaPRP,betaFR));
    end

    function beta = HyNG( self, g1, g0, d0 )
      % J. C. Gilbert and J. Nocedal
      % Global convergence properties of conjugate gradient methods for optimization,
      % SIAM J. Optim., 2 (1992), pp. 21–42.
      g0g0    = dot(g0,g0) ;
      betaPRP = dot(g1,g1-g0)/g0g0 ;
      betaFR  = dot(g1,g1)/g0g0 ;
      beta    = max(-betaFR,min(betaPRP,betaFR));
    end

    function beta = HyDY( self, g1, g0, d0 )
      % Y. H. Dai and Y. Yuan
      % An efficient hybrid conjugate gradient method for unconstrained optimization,
      % Ann. Oper. Res., 103 (2001), pp. 33–47.
      beta = dot(g1,g1)/max( dot(d0,g1-g0), -dot(g0,d0)) ;
    end

    function beta = betaEval( self, g1, g0, d0 )
      % Compute the value of `beta` for the selected nonlinear CG method
      switch ( self.method )
      case 'GRAD'  ; beta = 0 ; % standrd gradient method
      case 'HS'    ; beta = self.HS(g1,g0,d0) ;
      case 'HS+'   ; beta = self.HSplus(g1,g0,d0) ;
      case 'FR'    ; beta = self.FR(g1,g0,d0) ;
      case 'PRP'   ; beta = self.PRP(g1,g0,d0) ;
      case 'PRP+'  ; beta = self.PRPplus(g1,g0,d0) ;
      case 'CD'    ; beta = self.CD(g1,g0,d0) ;
      case 'LS'    ; beta = self.LS(g1,g0,d0) ;
      case 'DY'    ; beta = self.DY(g1,g0,d0) ;
      case 'N'     ; beta = self.N(g1,g0,d0) ;
      case 'HyTAS' ; beta = self.HyTAS(g1,g0,d0) ;
      case 'HyHS'  ; beta = self.HyHS(g1,g0,d0) ;
      case 'HyNG'  ; beta = self.HyNG(g1,g0,d0) ;
      case 'HyDY'  ; beta = self.HyDY(g1,g0,d0) ;
      otherwise
        error('MinimizationCG, betaEval no method selected (this should not be happen)');
      end
    end

  end

  methods

    function self = MinimizationCG( fun, ls )
      %
      % fun = function to be minimized
      % ls  = linesearch used, can be LinesearchArmijo, LinesearchWolfe, LinesearchGoldenSection
      self@MinimizationND( fun, ls ) ;
      self.method_names    = {'GRAD','HS','HS+','FR','PRP','PRP+','CD','LS','DY','N','HyTAS','HyHS','HyNG','HyDY'} ;
      self.method          = self.method_names{1} ;
      self.direction_short = 1e-3 ;
      self.angle_too_small = pi/180 ; % 1 degree
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
      % generic onjugate gradient minimization algorithm
      % x0 = initial point
      %
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
          d1   = -g1 ; % first iteration, search direction is -gradient
          beta = 0 ;
        else
          beta = self.betaEval( g1, g0, d0 ) ;
          d1   = -g1 + beta*d0 ;
          % check if the direction is descending, if direction is >= 89 degree from gradient reset to gradient directions
          % use >= to catch d1 == 0
          nrm_d1 = norm(d1) ;
          if dot( d1, g1 ) >= -nrm_d1 * nrm_g1 * (pi-self.angle_too_small)
            if self.debug_state
              fprintf(1,'direction angle about 90 degree, reset direction search, ...' ) ;
            end
            d1 = -g1 ; % reset direction
          elseif nrm_d1 <= self.direction_short * nrm_g1
            if self.debug_state
              fprintf(1,'direction length too short, reset direction search, ...' ) ;
            end
            d1 = -g1 ; % reset direction
          end
        end
      %
      if norm(d1,inf) == 0
        error('MinimizationCG, bad direction d == 0\n') ;
      end
        %
        % minimize along search direction
        [xs,alpha] = self.step1D( xs, d1, alpha ) ;
        %
        if self.debug_state
          fprintf(1,' alpha = %8.4g, beta = %8.4g\n', alpha, beta) ;
        end
        %
        % save old gradient and direction
        g0 = g1 ;
        d0 = d1 ;
      end
    end

  end
end
