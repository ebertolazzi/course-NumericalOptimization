%
% AA 2016/2017 Course: Numerical optimization
% by Enrico Bertolazzi
% email: enrico.bertolazzi@unitn.it
%
classdef LinesearchForwardBackward < handle
  % This class is the base class for a generic linesearch

  properties (SetAccess = protected, Hidden = true)
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
    barrier_reduce
    dumpMin      % minimum dumping factor
    dumpMax      % maximum dumping factor
    alpha_epsi   % minimum interval lenght for Wolfe linesearch
    debug_status % if true activate debug messages
    f0           % stored value f(0)
    Df0          % stored value f'(0)
    c1Df0        % stored value max(c1*f'(0),slopemax)
    name         % name of linesearch, set by the serived classed
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = LinesearchForwardBackward( name )
      %
      % constructor
      %
      self.c1             = 0.1;
      self.c2             = 0.2;
      self.tau_LS         = 1.1;
      self.tau_acc        = 1.2;
      self.alpha_min      = 1e-50;
      self.alpha_max      = 1e50;
      self.barrier_reduce = 1e-3;
      self.dumpMin        = 0.05;
      self.dumpMax        = 0.95;
      self.alpha_epsi     = eps^(1/3);
      self.debug_status   = false;
      self.name           = name;
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
    function y = ArmijoSlope( self, x )
      y = self.f0+self.c1Df0*x;
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
    %                     _          _   _
    %   __ _ _  _ __ _ __| |_ _ __ _| |_(_)__
    %  / _` | || / _` / _` | '_/ _` |  _| / _|
    %  \__, |\_,_\__,_\__,_|_| \__,_|\__|_\__|
    %     |_|
    %
    function alpha = quadratic( ~, f0, Df0, fp, p )
      alpha = Df0 * p^2 / ( 2*(f0+Df0*p-fp) );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %   ___                           _ ___          _                       _
    %  | __|__ _ ___ __ ____ _ _ _ __| | _ ) __ _ __| |____ __ ____ _ _ _ __| |
    %  | _/ _ \ '_\ V  V / _` | '_/ _` | _ \/ _` / _| / /\ V  V / _` | '_/ _` |
    %  |_|\___/_|  \_/\_/\__,_|_| \__,_|___/\__,_\__|_\_\ \_/\_/\__,_|_| \__,_|
    %
    function [LO,HI,ierr] = ForwardBackward( self, alpha_guess )
      % find alpha_min <= alpha0 < alpha1 <= alpha_max such that
      % alpha0 satify Armijo and alpha1 DO NOT satisfy Armijo
      % ierr =  0 interval found
      % ierr =  1 both alpha0, alpha1 satisfy Armijo f(alpha0) <  f(alpha1)
      % ierr =  2 both alpha0, alpha1 satisfy Armijo f(alpha0) >= f(alpha1)
      % ierr = -1 both alpha0, alpha1 DO NOTY satisfy Armijo
      % ierr = -2 Df0 >= 0
      % ierr = -3 f0 infinite or NaN
      % ierr = -4 step too small <= alpha_min

      sname = sprintf( 'Linesearch[%s]::ForwardBackward: ', self.name );

      % compute initial value and derivative
      f          = @(a) self.fun1D.eval(a);
      df         = @(a) self.fun1D.eval_D(a);
      self.f0    = f(0);
      self.Df0   = df(0);
      self.c1Df0 = 0;

      % check initial point
      if ~isfinite(self.f0)
        ierr     = -3;
        LO.alpha = 0;
        LO.f     = 0;
        LO.Df    = 0;
        HI       = LO;
        warning('%sf(0) = %g and f''(0) = %g must be a regular number\n', ...
                 sname, self.f0, self.Df0 );
        return;
      end

      % if not decreasing return an error
      if self.Df0 >= 0
        ierr     = -2;
        LO.alpha = 0;
        LO.f     = 0;
        LO.Df    = 0;
        HI       = LO;
        warning( '%sDf0 = %g must be negative\n', sname, LO.Df );
        return;
      end

      ierr = 0;

      % if fa1 DO NOT satisfy Armijo --> backward search
      self.c1Df0 = self.c1*self.Df0;
      HI.alpha   = 0;
      HI.f       = self.f0;
      HI.Df      = self.Df0;
      LO.alpha   = max(self.alpha_min,min(self.alpha_max,alpha_guess));
      LO.f       = f(LO.alpha);
      LO.Df      = df(LO.alpha);
      tauf       = self.tau_LS;
      armijo_ok  = false;
      while LO.alpha > self.alpha_min
        if isfinite(LO.f)
          armijo_ok = LO.f <= self.f0 + LO.alpha * self.c1Df0;
          if armijo_ok; break; end % satisfy Armijo condition
          if LO.f < self.f0 && LO.Df > 0
            break; % the interval contains a minimum
          end
        end
        % next point
        HI       = LO;
        LO.alpha = max(LO.alpha/tauf,self.alpha_min);
        LO.f     = f(LO.alpha);
        LO.Df    = df(LO.alpha);
        tauf     = tauf * self.tau_acc; % update tau factor
      end
      %
      if LO.alpha <= self.alpha_min
        ierr = -4;
        warning( '%sstep too small\n', sname );
        return;
      end
      %
      if armijo_ok
        if HI.alpha < LO.alpha
          % satisfy armijo at first step, try to enlarge interval
          HI       = LO;
          LO.alpha = 0;
          LO.f     = self.f0;
          LO.Df    = self.Df0;
          tauf     = self.tau_LS;
          %
          N = HI;
          while N.alpha < self.alpha_max && ...
                (N.f <= self.f0 + N.alpha * self.c1Df0) && ...
                (N.Df+self.c1Df0 < 0)
            LO       = HI;
            HI       = N;
            N.alpha  = tauf * N.alpha;
            N.f      = f(N.alpha);
            N.Df     = df(N.alpha);
            tauf     = tauf * self.tau_acc; % update tau factor
          end
          %
          return;
        elseif LO.Df >= 0
          % minimum in [0,LO.alpha]
          HI       = LO;
          LO.alpha = 0;
          LO.f     = self.f0;
          LO.Df    = self.Df0;
          return;
        elseif ~isfinite(HI.f)
          % minimum in [LO.alpha,HI.alpha] but interval must be reduced
          dAlphaMin = max( self.barrier_reduce * (HI.alpha-LO.alpha), self.alpha_min );
          % search for a finite value
          while ~isfinite(HI.f) && HI.alpha-LO.alpha > dAlphaMin
            tmp.alpha = (LO.alpha+HI.alpha)/2;
            tmp.f     = f(tmp.alpha);
            tmp.Df    = df(tmp.alpha);
            if isfinite(tmp.f) && ...
               tmp.f <= self.f0 + tmp.alpha * self.c1Df0 && ...
               tmp.Df < 0
              LO = tmp;
            else
              HI = tmp;
            end
          end
        else
          % HI.f is finite and LO.Df < 0 minimum in [LO.alpha,HI.alpha]
        end
        return;
      end
      %
      % at this point
      % LO.f < self.f0 && LO.Df > 0;
      % minimum in [0,LO.alpha]
      HI       = LO;
      LO.alpha = 0;
      LO.f     = self.f0;
      LO.Df    = self.Df0;

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
