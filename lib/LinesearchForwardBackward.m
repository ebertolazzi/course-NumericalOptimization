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

    n_fun_eval
    max_fun_eval % a positive integer input variable.
                 % Linesearch termination occurs when the number of calls
                 % to f(alpha) is at least max_fun_eval by the end of an iteration.

    c1           % Armijo constant to accept the step (0,1/2]: f(alpha) <= f(0) + c1 * alpha * f'(0)
    c2           % Wolfe constant to accept the step [c1,1/2]: f'(alpha) >= c2 * f'(0)

    tau_LS       % multiplicative factor for Forward search
    tau_acc      % modify tau to accelerate exploration for large or very small interval

    alpha_min    % minimum accepted step
    alpha_max    % maximum accepted step

    barrier_reduce
    debug_status % if true activate debug messages

    f0           % stored value f(0)
    Df0          % stored value f'(0)
    c1Df0        % stored value max(c1*f'(0),slopemax)
    c2AbsDf0     % stored value c2*abs(f'(0))

    name         % name of linesearch, set by the serived classed
  end

  methods (Hidden = true)
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [L,R] = infGrow( self, L, R, alpha_guess )
      % minimum in [0,alpha_guess] but interval must be reduced
      dAlphaMin = max( self.barrier_reduce * alpha_guess, self.alpha_min );
      % search for a finite value
      L0 = L;
      while ~isfinite(R.f)
        tmp.alpha         = (L.alpha+R.alpha)/2;
        [ tmp.f, tmp.Df ] = self.fDf( tmp.alpha );
        if isfinite(tmp.f) && ...
           tmp.f <= self.f0 + tmp.alpha * self.c1Df0 && ...
           tmp.Df < 0
          L = tmp;
        else
          R = tmp;
        end
        if R.alpha-L.alpha < dAlphaMin && L.alpha > 0; break; end
      end

      if ~isfinite(R.f)
        R = L;
        L = L0;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,Df] = fDf( self, alpha )
      self.n_fun_eval = self.n_fun_eval+1;
      f  = self.fun1D.eval( alpha );
      Df = self.fun1D.eval_D( alpha );
    end
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
      self.max_fun_eval   = 1000;
      self.debug_status   = false;
      self.name           = name;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setFunction( self, fun1D )
      % set the function object used in the 1D minimization
      self.fun1D = fun1D;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setInitialTau( self, tau_LS )
      % set the initial dumping factor tau used in forward backward search
      if ( tau_LS > 1 ) || (tau_LS < 1000 )
        error('Linesearch[%s]::setInitialTau, constant tau_LS = %g must be > 1 and < 1000, tau_LS = %g\n',self.name,tau_LS);
      end
      self.tau_LS = tau_LS;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setAccelerationTau( self, tau_acc )
      % set acceleration factor used in forward backward search
      if ( tau_acc >= 1 ) || (tau_acc < 10 )
        error('Linesearch[%s]::setAccelerationTau, acceleration factor tau = %g must be >= 1 and < 10, tau_acc = %g\n',self.name,tau_acc);
      end
      self.tau_acc = tau_acc;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setC1C2( self, c1, c2 )
      % set the coefficients c1 (Armijo) and c2 (Wolfe) for lineasearch
      seps = sqrt(eps);
      % set the c1 coefficients for lineasearch
      if c1 > 0.5
        fprintf('Linesearch[%s], constant c1 = %g must be <= 0.5, set to 0.5\n',self.name,c1);
        self.c1 = 0.5;
      elseif c1 < seps
        fprintf('Linesearch[%s], constant c1 = %g must be >= %g, set to %g\n',self.name,c1,seps,seps);
        self.c1 = seps;
      else
        self.c1 = c1;
      end
      % set the c2 coefficients for lineasearch
      if c2 < self.c1
        fprintf('Linesearch[%s], constant c2 = %g must be >= c1 = %g\n',self.name,c2,self.c1);
        self.c2 = self.c1;
      elseif c2 > 0.5
        fprintf('Linesearch[%s], constant c2 = %g must be <= 0.5, set to 0.5\n',self.name,c2);
        self.c2 = 0.5;
      else
        self.c2 = c2;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = ArmijoSlope( self, x )
      y = self.f0+self.c1Df0*x;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function debug_on( self )
      self.debug_status = true;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function debug_off( self )
      self.debug_status = false;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotDebug( self, alpha_guess )
      figure();
      subplot(3,1,1);
      self.plot(alpha_guess*10);
      subplot(3,1,2);
      self.plot(alpha_guess);
      subplot(3,1,3);
      self.plot(alpha_guess/10);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %                     _          _   _
    %   __ _ _  _ __ _ __| |_ _ __ _| |_(_)__
    %  / _` | || / _` / _` | '_/ _` |  _| / _|
    %  \__, |\_,_\__,_\__,_|_| \__,_|\__|_\__|
    %     |_|
    %
    function alpha = quadratic( ~, f0, Df0, fp, p )
      alpha = Df0 * p^2 / ( 2*(f0+Df0*p-fp) );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %   ___                           _ ___          _                       _
    %  | __|__ _ ___ __ ____ _ _ _ __| | _ ) __ _ __| |____ __ ____ _ _ _ __| |
    %  | _/ _ \ '_\ V  V / _` | '_/ _` | _ \/ _` / _| / /\ V  V / _` | '_/ _` |
    %  |_|\___/_|  \_/\_/\__,_|_| \__,_|___/\__,_\__|_\_\ \_/\_/\__,_|_| \__,_|
    %
    function [LO,HI,ierr] = ForwardBackward( self, alpha_guess )
      % find alpha_min <= alpha0 < alpha1 <= alpha_max such that
      % alpha0 satify Armijo and alpha1 DO NOT satisfy Armijo
      % ierr =  0 interval found
      %
      % ierr = -1 f(0) infinite or NaN
      % ierr = -2 f'(0) >= 0
      self.n_fun_eval = 0;
      LO = {};
      HI = {};
      sname = sprintf( 'Linesearch[%s]::ForwardBackward: ', self.name );

      %
      % compute initial value and derivative
      %
      [ self.f0, self.Df0 ] = self.fDf(0);
      self.c1Df0    = 0;
      self.c2AbsDf0 = 0;

      %
      % check initial point
      if ~isfinite(self.f0)
        ierr = -1;
        fprintf('%sf(0) = %g and f''(0) = %g must be a regular number\n', ...
                 sname, self.f0, self.Df0 );
        return;
      end

      %
      % if not decreasing return an error
      if self.Df0 >= 0
        ierr = -2;
        fprintf( '%sDf0 = %g must be negative\n', sname, self.Df0 );
        return;
      end

      %
      % in case f(alpha) is infinite try to detect the barrier
      L.alpha       = 0;
      L.f           = self.f0;
      L.Df          = self.Df0;
      R.alpha       = alpha_guess;
      [ R.f, R.Df ] = self.fDf(alpha_guess);
      %
      % check if f(alpha) is infinite
      if ~isfinite(R.f)
        [L,R] = self.infGrow( L, R, alpha_guess );
      end

      %
      % if step too small exit
      if R.alpha <= self.alpha_min
        ierr = -3;
        fprintf(2,'\n%sstep too small\n', sname );
        return;
      end

      %
      % check Armijo
      self.c1Df0    = self.c1*self.Df0;
      self.c2AbsDf0 = self.c2*abs(self.Df0);
      if (R.f - self.f0) > R.alpha*self.c1Df0
        % Armijo NOT satified
        % reduce the step until f(R) <= f(L) < f(0)
        tauf          = self.tau_LS;
        L.alpha       = max(R.alpha/tauf,self.alpha_min);
        [ L.f, L.Df ] = self.fDf( L.alpha );
        S             = R;
        Llesscnt      = 0;
        while R.alpha > self.alpha_min
          if L.f < self.f0
            if L.f >= R.f; break; end
            Llesscnt = Llesscnt+1;
            if Llesscnt > 10; break; end
          end
          S             = R;
          R             = L;
          L.alpha       = max(R.alpha/tauf,self.alpha_min);
          [ L.f, L.Df ] = self.fDf( L.alpha );
          tauf          = tauf * self.tau_acc; % update tau factor
        end
        if L.f < self.f0
          R = S;
        else
          ierr = -3;
          fprintf(2,'\n%sstep too small 2\n', sname );
          return;
        end
      else
        % satisfy Armijo at first step, try to enlarge interval
        tauf          = self.tau_LS;
        N.alpha       = tauf * R.alpha;
        [ N.f, N.Df ] = self.fDf( N.alpha );
        while N.alpha < self.alpha_max && ...
              (N.f-self.f0) <= N.alpha * self.c1Df0 && ...
              N.f <= R.f && N.Df < 0
          L             = R;
          R             = N;
          N.alpha       = tauf * N.alpha;
          [ N.f, N.Df ] = self.fDf( N.alpha );
          tauf          = tauf * self.tau_acc; % update tau factor
        end

        if ~isfinite(N.f)
          [L,R] = self.infGrow( L, N, alpha_guess );
        end
      end
      %
      % at this point
      % L.f < self.f0 && L.Df > 0;
      % minimum in [0,L.alpha]
      ierr    = 0;
      if L.f > R.f
        LO = R;
        HI = L;
      else
        LO = L;
        HI = R;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plot( self, alpha_max )
      alpha = -alpha_max/10:alpha_max/1000:alpha_max;
      y     = zeros(size(alpha));
      for k=1:length(alpha)
        y(k) = self.fun1D.eval(alpha(k));
      end
      plot( alpha, y );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
