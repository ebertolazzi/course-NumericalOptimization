classdef LinesearchForwardBackward < handle
  % This class is the base class for a generic linesearch

  properties (SetAccess = private, Hidden = true)
    fun1D   % Object used in the linesearch that contains the methods
            % f(alpha)   = fun1D.eval(alpha) 
            % f'(alpha)  = fun1D.eval_D(alpha) 
            % f''(alpha) = fun1D.eval_D_DD(alpha)
    c1      % Armijo constant to accept the step (0,1/2]
            % f(alpha) <= f(0) + c1 * alpha * f'(0)
    tau_LS  % multiplicative factor for Forward search
    tau_acc % modify tau to accelerate exploration for large or very small interval
    alpha_min
    alpha_max
    debug_on
  end
  
  methods

    function self = LinesearchForwardBackward()
      self.c1        = 0.1 ;
      self.tau_LS    = 1.5 ;
      self.tau_acc   = 1.2 ;
      self.alpha_min = 1e-10 ;
      self.alpha_max = 1e10 ;
      self.debug_on  = true ;
    end

    function self = setFunction( self, fun1D )
      % set the function object used in the 1D minimization
      self.fun1D = fun1D ;
    end

    function self = setInitialTau( self, tau_LS )
      % set the initial dumping factor tau
      if ( tau_LS > 1 ) || (tau_LS < 1000 )
        error('LinesearchForwardBackward::setInitialTau, constant tau_LS = %g must be > 1 and < 1000, tau_LS = %g\n',tau_LS);
      end
      self.tau_LS = tau_LS ;
    end

    function self = setAccelerationTau( self, tau_acc )
      % set the initial dumping factor tau
      if ( tau_acc > 1 ) || (tau_acc < 10 )
        error('LinesearchForwardBackward::setAccelerationTau, acceleration factor tau = %g must be > 1 and < 10, tau_acc = %g\n',tau_acc);
      end
      self.tau_acc = tau_acc ;
    end

    function self = setAlphaRange( self, arange )
      % set the range for alpha step
      if length(arange) ~= 2
        error('LinesearchForwardBackward::setAlphaRange, argument must be a 2 element vector\n');
      end
      if arange(1) >= arange(2) && arange(1) > 1e-50 && arange(2) < 1e50
        error('LinesearchForwardBackward::setAlphaRange, bad range [%g,%g]\n',arange(1),arange(2));
      end
      self.alpha_min = arange(1);
      self.alpha_max = arange(2);
    end

    function self = setC1( self, c1 )
      % set the c1 coefficients for lineasearch
      if (c1 >= 0.5) || (c1 < sqrt(eps) )
        error('LinesearchForwardBackward, constant c1 = %g must be <= 0.5 and >= %10f\n',c1,sqrt(eps));
      end
      self.c1 = c1 ;
    end

    function self = setDebug( self )
      self.debug_on = true ;
    end

    function self = setNoDebug( self )
      self.debug_on = false ;
    end

    function [alpha0,alpha1,ierr] = ForwardBackward( self, alpha_guess )
      % find alpha_min <= alpha0 < alpha1 <= alpha_max such that 
      % alpha0 satify Armijo and alpha1 DO NOT satisfy Armijo
      % ierr = 0  interval found
      % ierr = 1  both alpha0, alpha1 satisfy Armijo
      % ierr = -1 both alpha0, alpha1 DO NOTY satisfy Armijo
      % ierr = -2 Df0 >= 0
      if (alpha_guess <= 0) || (alpha_guess >= self.alpha_max)
        error('LinesearchForwardBackward::ForwardBackward, bad alpha_guess = %g with alpha_max = %g\n', alpha_guess, self.alpha_max );        
      end
      Df0 = self.fun1D.eval_D(0) ;
      if (Df0 >= 0)
        figure();
        subplot(3,1,1);
        self.plot(alpha_guess*10);
        subplot(3,1,2);
        self.plot(alpha_guess);
        subplot(3,1,3);
        self.plot(alpha_guess/10);
        error('LinesearchForwardBackward::ForwardBackward, Df0 = %g must be negative\n', Df0 );        
      end
      f0   = self.fun1D.eval(0) ;
      fa   = self.fun1D.eval(alpha_guess) ;
      tauf = self.tau_LS ;
      ierr = 0 ;
      if fa <= f0 + alpha_guess * self.c1 * Df0
        % satisfy Armijo --> forward search
        alpha0 = alpha_guess ;
        alpha1 = tauf*alpha0 ;
        while self.fun1D.eval(alpha1) <= f0 + alpha1 * self.c1 * Df0
          alpha0 = alpha1 ;
          alpha1 = tauf*alpha0 ;
          if alpha1 > self.alpha_max ; ierr = 1 ; break ; end
          tauf = tauf * self.tau_acc ; % update tau factor
        end
        if alpha1 > self.alpha_max
          alpha1 = self.alpha_max ;
          if self.fun1D.eval(alpha1) > f0 + alpha1 * self.c1 * Df0
            ierr = 0; % fond             
          end          
        end        
      else
        % DO NOT satisfy Armijo --> backward search
        alpha1 = alpha_guess ;
        alpha0 = alpha1/tauf ;
        while self.fun1D.eval(alpha0) > f0 + alpha1 * self.c1 * Df0
          alpha1 = alpha0 ;
          alpha0 = alpha1/tauf ;
          if alpha0 < self.alpha_min ; ierr = -1 ; break ; end
          tauf = tauf * self.tau_acc ; % update tau factor
        end
        if alpha0 < self.alpha_min
          alpha0 = self.alpha_min ;
          if self.fun1D.eval(alpha0) <= f0 + alpha0 * self.c1 * Df0
            ierr = 0; % found             
          end          
        end
      end
    end
    
    function plot( self, alpha_max )
      alpha = -alpha_max/10:alpha_max/1000:alpha_max ;
      y = [] ;
      for a=alpha ; y = [y self.fun1D.eval(a)]; end
      plot( alpha, y );
    end

  end
end
