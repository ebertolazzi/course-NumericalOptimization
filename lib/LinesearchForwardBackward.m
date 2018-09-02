%
% AA 2016/2017 Course: Numerical optimization
% by Enrico Bertolazzi
% email: enrico.bertolazzi@unitn.it
%
classdef LinesearchForwardBackward < handle
  % This class is the base class for a generic linesearch

  properties (SetAccess = private, Hidden = true)
    fun1D        % Object used in the linesearch that contains the methods
                 % f(alpha)   = fun1D.eval(alpha)
                 % f'(alpha)  = fun1D.eval_D(alpha)
                 % f''(alpha) = fun1D.eval_D_DD(alpha)
    c1           % Armijo constant to accept the step (0,1/2]: f(alpha) <= f(0) + c1 * alpha * f'(0)
    c2           % Wolfe constant to accept the step [c1,1/2]: f'(alpha) >= c2 * f'(0)
    tau_LS       % multiplicative factor for Forward search
    tau_acc      % modify tau to accelerate exploration for large or very small interval
    alpha_min    % minimum accepted step
    alpha_max    % maximum accepted step
    dumpMin      % minimum dumping factor
    dumpMax      % maximum dumping factor
    alpha_epsi   % minimum interval lenght for Wolfe linesearch
    debug_status % if true activate debug messages
    f0           % stored value f(0)
    Df0          % stored value f'(0)
    name         % name of linesearch, set by the serived classed
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = LinesearchForwardBackward( name )
      %
      % constructor
      %
      self.c1           = 0.05;
      self.c2           = 0.2;
      self.tau_LS       = 1.5;
      self.tau_acc      = 1.2;
      self.alpha_min    = 1e-50;
      self.alpha_max    = 1e50;
      self.dumpMin      = 0.05;
      self.dumpMax      = 0.95;
      self.alpha_epsi   = eps^(1/3);
      self.debug_status = false;
      self.name         = name;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setFunction( self, fun1D )
      % set the function object used in the 1D minimization
      self.fun1D = fun1D;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setInitialTau( self, tau_LS )
      % set the initial dumping factor tau used in forward backward search
      if ( tau_LS > 1 ) || (tau_LS < 1000 )
        error('Linesearch[%s]::setInitialTau, constant tau_LS = %g must be > 1 and < 1000, tau_LS = %g\n',self.name,tau_LS);
      end
      self.tau_LS = tau_LS;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setAccelerationTau( self, tau_acc )
      % set acceleration factor used in forward backward search
      if ( tau_acc >= 1 ) || (tau_acc < 10 )
        error('Linesearch[%s]::setAccelerationTau, acceleration factor tau = %g must be >= 1 and < 10, tau_acc = %g\n',self.name,tau_acc);
      end
      self.tau_acc = tau_acc;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setAlphaRange( self, amin, amax )
      % set the range for alpha step
      if ~ (isscalar(amin) && isscalar(amax))
        error('Linesearch[%s]::setAlphaRange, arguments must be a 2 scalars\n',self.name);
      end
      if amin >= amax && amin > 1e-50 && amax < 1e50
        error('Linesearch[%s]::setAlphaRange, bad range [%g,%g]\n',self.name,amin,amax);
      end
      self.alpha_min = amin;
      self.alpha_max = amax;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setC1C2( self, c1, c2 )
      % set the coefficients c1 (Armijo) and c2 (Wolfe) for lineasearch
      seps = sqrt(eps);
      % set the c1 coefficients for lineasearch
      if c1 > 0.5
        warning('Linesearch[%s], constant c1 = %g must be <= 0.5, set to 0.5\n',self.name,c1);
        self.c1 = 0.5;
      elseif c1 < seps
        warning('Linesearch[%s], constant c1 = %g must be >= %g, set to %g\n',self.name,c1,seps,seps);
        self.c1 = seps;
      else
        self.c1 = c1;
      end
      % set the c2 coefficients for lineasearch
      if c2 < self.c1
        warning('Linesearch[%s], constant c2 = %g must be >= c1 = %g\n',self.name,c2,self.c1);
        self.c2 = self.c1;
      elseif c2 > 0.5
        warning('Linesearch[%s], constant c2 = %g must be <= 0.5, set to 0.5\n',self.name,c2);
        self.c2 = 0.5;
      else
        self.c2 = c2;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setDump( self, d1, d2 )
      % set dumping coefficients for Zoom
      if d1 > 0.4
        warning('Linesearch[%s], dump d1 = %g must be <= 0.4, set to 0.4\n',self.name,d1);
        self.dumpMin = 0.4;
      else
        self.dumpMin = d1;
      end
      if d2 < 0.6
        warning('Linesearch[%s], dump d2 = %g must be >= 0.6, set to 0.6\n',self.name,d2);
        self.dumpMax = 0.6;
      else
        self.dumpMax = d2;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function debug_on( self )
      self.debug_status = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function debug_off( self )
      self.debug_status = false;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotDebug( self, alpha_guess )
      figure();
      subplot(3,1,1);
      self.plot(alpha_guess*10);
      subplot(3,1,2);
      self.plot(alpha_guess);
      subplot(3,1,3);
      self.plot(alpha_guess/10);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %                         _           _   _
    %    __ _ _   _  __ _  __| |_ __ __ _| |_(_) ___
    %   / _` | | | |/ _` |/ _` | '__/ _` | __| |/ __|
    %  | (_| | |_| | (_| | (_| | | | (_| | |_| | (__
    %   \__, |\__,_|\__,_|\__,_|_|  \__,_|\__|_|\___|
    %      |_|
    %
    function alpha = quadratic( ~, f0, Df0, fp, p )
      alpha = Df0 * p^2 / ( 2*(f0+Df0*p-fp) );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [alpha0,alpha1,fa0,fa1,ierr] = ForwardBackward( self, alpha_guess )
      % find alpha_min <= alpha0 < alpha1 <= alpha_max such that
      % alpha0 satify Armijo and alpha1 DO NOT satisfy Armijo
      % ierr = 0  interval found
      % ierr = 1  both alpha0, alpha1 satisfy Armijo f(alpha0) <  f(alpha1)
      % ierr = 2  both alpha0, alpha1 satisfy Armijo f(alpha0) >= f(alpha1)
      % ierr = -1 both alpha0, alpha1 DO NOTY satisfy Armijo
      % ierr = -2 Df0 >= 0

      % correct alpha_guess into the required interval
      alpha0   = max(self.alpha_min,min(self.alpha_max,alpha_guess));
      % compute initial value and derivative
      f        = @(a) self.fun1D.eval(a);
      df       = @(a) self.fun1D.eval_D(a);
      self.f0  = f(0);
      self.Df0 = df(0);
      % if not decreasing issue an error (if in debug also plot the function)
      if (self.Df0 >= 0)
        if self.debug_status
          self.plotDebug(alpha_guess);
          error('Linesearch[%s]::ForwardBackward, Df0 = %g must be negative\n',self.name,self.Df0 );
        else
          ierr   = -2;
          alpha0 = 0;
          alpha1 = 0;
          fa0    = 0;
          fa1    = 0;
          warning('Linesearch[%s]::ForwardBackward, Df0 = %g must be negative\n',self.name,self.Df0 );
        end
        return;
      end
      % if Df0 too big cut it
      self.Df0 = max( self.Df0, -1e10 );

      % initialize search parameters
      c1Df0 = self.c1*self.Df0;
      tauf  = self.tau_LS;
      ierr  = 0;
      % decide if do forward or backward search
      fa0 = f(alpha0);
      if (fa0-self.f0) <= alpha0 * c1Df0
        % satisfy Armijo --> forward search

        % if increasing minima in [0,alpha0]
        if df(alpha0) >= 0
          alpha1 = alpha0; fa1 = fa0;
          alpha0 = 0;      fa0 = self.f0;
          ierr = 1;
          return;
        end

        % not increasing search next interval
        alpha1 = tauf * alpha0; fa1 = f(alpha1);
        % continue to loop until satify Armijo and is non increasing
        while true
          % check if found point that violates Armijo
          if (fa1-self.f0) > alpha1 * c1Df0; return; end
          % if increasing break
          if fa1 > fa0;       ierr = 1; return; end
          if df(alpha1) >= 0; ierr = 2; return; end

          % check if interval become too large
          if tauf * alpha0 >= self.alpha_max
            % last interval to check
            alpha0 = alpha1;         fa0 = fa1;
            alpha1 = self.alpha_max; fa1 = f(alpha1);
            if fa1 > self.f0 + alpha1 * c1Df0; return; end % found!
            if fa1 > fa0; ierr = 1; else ierr = 2; end
            return;
          end
          % prepare for next interval
          alpha0 = alpha1;        fa0 = fa1;
          alpha1 = tauf * alpha0; fa1 = f(alpha1);
          tauf   = tauf * self.tau_acc; % update tau factor
        end
        % never pass here
      else
        % DO NOT satisfy Armijo --> backward search
        alpha1 = alpha_guess; fa1 = fa0;
        alpha0 = alpha1/tauf; fa0 = f(alpha0);
        while (fa0-self.f0) > alpha0 * c1Df0 % DO NOT satisfy Armijo
          alpha1 = alpha0; fa1 = fa0;
          if alpha1 <= self.alpha_min * tauf
            % force termination, last check
            alpha0 = self.alpha_min; fa0 = f(alpha0);
            if (fa0-self.f0) > alpha0 * c1Df0; ierr = -1; end
            return;
          end
          alpha0 = alpha1/tauf; fa0 = f(alpha0);
          tauf   = tauf * self.tau_acc; % update tau factor
        end
        % exiting from loop satisfying Armijo
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [alpha,ok] = Zoom( self, aLO, fLO, aHI, fHI, strongWolfe )
      % refine approximate minimum until it satisfy Wolfe condition

      f  = @(a) self.fun1D.eval(a);
      df = @(a) self.fun1D.eval_D(a);

      self.f0  = self.fun1D.eval(0.0);
      self.Df0 = self.fun1D.eval_D(0.0);

      c1Df0 = self.c1 * self.Df0;
      c2Df0 = self.c2 * self.Df0;

      % check if aLo satisfy wolfe condition
      DfLO = df(aLO);
      if DfLO >= c2Df0 && ( ~strongWolfe || DfLO <= -c2Df0 )
        alpha = aLO;
        ok    = true;
        return;
      end

      %   ____  _____ _____ ___ _   _ _____
      %  |  _ \| ____|  ___|_ _| \ | | ____|
      %  | |_) |  _| | |_   | ||  \| |  _|
      %  |  _ <| |___|  _|  | || |\  | |___
      %  |_| \_\_____|_|   |___|_| \_|_____|
      %

      Delta    = aHI - aLO;
      minDelta = Delta * self.alpha_epsi;

      while abs(Delta) > minDelta
        deltaLambda = self.quadratic( fLO, DfLO, fHI, Delta );
        if Delta > 0
          deltaLambda = min( Delta*self.dumpMax, max( Delta*self.dumpMin, deltaLambda ));
        else
          deltaLambda = max( Delta*self.dumpMax, min( Delta*self.dumpMin, deltaLambda ));
        end

        alpha = aLO + deltaLambda; fa = f(alpha);

        if fa > self.f0 + alpha*c1Df0 || fa > fLO
          aHI = alpha; fHI = fa;
        else
          Dfa = df(alpha);
          if Dfa >= c2Df0 && ( ~strongWolfe || Dfa <= -c2Df0 )
            ok = true;
            return; % found Wolfe point
          end
          % choose left or right interval
          if Dfa*(aLO-alpha) < 0
            aHI = aLO; fHI = fLO;
          end
          aLO = alpha; fLO = fa; DfLO = Dfa;
        end
        Delta = aHI - aLO;
      end
      ok    = false;
      alpha = aLO;
      warning( 'Linesearch[%s] (Zoom): failed aHI=%g aLO=%g DfLO=%g\n', self.name, aHI, aLO, DfLO );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, alpha_max )
      alpha = -alpha_max/10:alpha_max/1000:alpha_max;
      y = [];
      for a=alpha; y = [y self.fun1D.eval(a)]; end
      plot( alpha, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
