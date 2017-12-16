classdef LinesearchArmijo < handle
  % This class is an abstract class for a generic linesearch

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
  end
  
  methods

    function self = LinesearchArmijo()
      self.c1        = 0.1 ;
      self.tau_LS    = 1.5 ;
      self.tau_acc   = 1.2 ;
      self.alpha_min = 1e-10 ;
      self.alpha_max = 1e10 ;
    end

    function self = setC1( self, c1 )
      if (c1 >= 0.5) || (c1 < sqrt(eps) )
        error('Linesearch1D, constant c1 = %g must be <= 0.5 and >= %10f\n',c1,sqrt(eps));
      end
      self.c1 = c1 ;
    end
    
    % aggiungere altri parametri da cambiare

    function self = setFunction( self, fun1D )
      self.fun1D = fun1D ;
    end

    function [alpha_ott,ok] = search( self, alpha_guess )
      [alpha0,~,ierr] = self.ForwardBackward( alpha_guess );
      alpha_ott = alpha0 ;
      ok        = ierr >= 0 ;
    end

    function [alpha0,alpha1,ierr] = ForwardBackward( self, alpha_guess )
      % find alpha_min <= alpha0 < alpha1 <= alpha_max such that 
      % alpha0 satify Armijo and alpha1 DO NOT satisfy Armijo
      % ierr = 0  interval found
      % ierr = 1  both alpha0, alpha1 satisfy Armijo
      % ierr = -1 both alpha0, alpha1 DO NOTY satisfy Armijo
      % ierr = -2 Df0 >= 0
      if (alpha_guess <= 0) || (alpha_guess >= self.alpha_max)
        error('Linesearch1D::ForwardBackward, bad alpha_guess = %g with alpha_max = %g\n', alpha_guess, self.alpha_max );        
      end      
      Df0 = self.fun1D.eval_D(0) ;
      if (Df0 >= 0)
        error('Linesearch1D::ForwardBackward, Df0 = %g must be negative\n', Df0 );        
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
            ierr = 0; % fond             
          end          
        end
      end
    end
    
    function plot( self, alpha_max )
      alpha = 0:alpha_max/1000:alpha_max ;
      y = [] ;
      for a=alpha ; y = [y self.fun1D.eval(a)]; end
      plot( alpha, y );
    end
  end
end
