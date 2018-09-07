classdef LinesearchMoreThuente < LinesearchForwardBackward
  %%
  %  Translation of minpack subroutine cvsrch
  %  based on matlab conversiuon by Dianne O'Leary July 1991
  %  original code by
  %
  %  by Jorge J. More' and David J. Thuente
  %  Argonne National Laboratory. MINPACK Project. June 1983
  %
  %  OO matlab version by Enrico Bertolazzi
  %
  %  The purpose of the linesearch is to find a step which satisfies
  %  a sufficient decrease condition and a curvature condition.
  %  The user must provide a subroutine which calculates the
  %  function and the gradient.
  %
  %  At each stage the method `search` updates an interval of
  %  uncertainty with enM.Dfoints LO.alpha and HI.alpha.
  %  The interval of uncertainty is initially chosen so that it
  %  contains a minimizer of the modified function
  %
  %    phi(alpha) = f(alpha) - f(0) - c1*alpha*f'(0)
  %
  %  If a step is obtained for which the modified function has a nonpositive
  %  function value and nonnegative derivative, then the interval of
  %  uncertainty is chosen so that it contains a minimizer of f(alpha).
  %
  %  The algorithm is designed to find a step which satisfies
  %  the sufficient decrease condition
  %
  %    f(alpha) <= f(0) + c1*alpha*f'(0),
  %
  %  and the curvature condition
  %
  %    | f'(alpha) | <= c2 * |f'(0)|
  %
  %  If c1 is less than c2 and if, for example, the function is bounded below,
  %  then there is always a step which satisfies both conditions.
  %  If no step can be found which satisfies both conditions, then the
  %  algorithm usually stops when rounding errors prevent further progress.
  %  In this case alpha only satisfies the sufficient decrease condition.
  %

  properties (SetAccess = private, Hidden = true)
    xtol         % Linesearch termination occurs when the relative width
                 % of the interval of uncertainty is at most xtol.
    max_fun_eval % a positive integer input variable.
                 % Linesearch termination occurs when the number of calls
                 % to f(alpha) is at least max_fun_eval by the end of an iteration.
    info         % an integer output variable set as follows:
                 % info = 0  Improper input parameters.
                 %
                 % info = 1  The sufficient decrease condition and the
                 %           directional derivative condition hold.
                 %
                 % info = 2  Relative width of the interval of uncertainty
                 %           is at most xtol.
                 %
                 % info = 3  Number of calls to fcn has reached max_fun_eval.
                 %
                 % info = 4  The step is at the lower bound alpha_min.
                 %
                 % info = 5  The step is at the upper bound alpha_max.
                 %
                 % info = 6  Rounding errors prevent further progress.
                 %           There may not be a step which satisfies the
                 %           sufficient decrease and curvature conditions.
                 %           Tolerances may be too small.
    xtrapf
  end

  methods (Hidden = true)

    function [LO, HI, alphaf, bracketed, info] = safeguardedStep( self, LO, HI, M, bracketed )
      %
      % The purpose of cstep is to compute a safeguarded step for
      % a linesearch and to update an interval of uncertainty for
      % a minimizer of the function.
      %
      % The struct LO contains the step with the least function value.
      % The struct M contains the current step.
      % It is assumed that the derivative at LO is negative in the direction
      % of the step.
      % If bracketed is set true then a minimizer has been bracketed in an
      % interval of uncertainty with endpoints LO.alpha and HI.alpha.
      %
      % LO specify the step at the best step obtained so far.
      % The derivative must be negative in the direction of the step, that is,
      % LO.Df and (M.alpha-LO.alpha) must have opposite signs.
      % On output these parameters are updated appropriately.
      %
      % HI specify the step at the other endpoint of the interval of uncertainty.
      % output these parameters are updated appropriately.
      %
      % struct M specify the current step.
      % If bracketed is set true then on input M.alpha must be between
      % LO.alpha and HI.alpha. On output M.alpha is set to the new step.
      %
      % bracketed is a logical variable which specifies if a minimizer
      % has been bracketed. If the minimizer has not been bracketed
      % then on input bracketed must be set false.
      % If the minimizer is bracketed then on output bracketed is set true.
      %
      % info is an integer output variable set as follows:
      % If info = 1,2,3,4,5, then the step has been computed according to one
      % of the five cases below.
      % Otherwise info = 0, and this indicates improper input parameters.
      %
      % Check the input parameters for errors.
      aL = min(LO.alpha,HI.alpha);
      aR = max(LO.alpha,HI.alpha);
      if bracketed
        if M.alpha <= aL || M.alpha >= aR
          info = 0;
          return;
        end
      end
      if LO.Df*(M.alpha-LO.alpha) >= 0
        info = 0;
        return;
      end
      %
      % Determine if the derivatives have opposite sign.
      %
      sgnd = M.Df*sign(LO.Df);
      %
      % First case.
      % A higher function value. The minimum is bracketed.
      % If the cubic step is closer to LO.alpha than the quadratic step,
      % the cubic step is taken, else the average of the cubic and
      % quadratic steps is taken.
      %
      if M.f > LO.f
        info  = 1;
        bound = true;
        theta = 3*(LO.f - M.f)/(M.alpha - LO.alpha) + LO.Df + M.Df;
        s     = max(abs([theta,LO.Df,M.Df]));
        gamma = s*sqrt((theta/s)^2 - (LO.Df/s)*(M.Df/s));
        if M.alpha < LO.alpha
          gamma = -gamma;
        end
        p = (gamma - LO.Df) + theta;
        q = ((gamma - LO.Df) + gamma) + M.Df;
        r = p/q;
        alphac = LO.alpha + r*(M.alpha - LO.alpha);
        alphaq = LO.alpha + ((LO.Df/((LO.f-M.f)/(M.alpha-LO.alpha)+LO.Df))/2)*(M.alpha - LO.alpha);
        if abs(M.alphac-LO.alpha) < abs(M.alphaq-LO.alpha)
          alphaf = alphac
        else
          alphaf = alphac + (alphaq - alphac)/2;
        end
        bracketed = true;
      elseif sgnd < 0
        %
        % Second case.
        % A lower function value and derivatives of opposite sign.
        % The minimum is bracketed. If the cubic step is closer to LO.alpha
        % than the quadratic (secant) step, the cubic step is taken,
        % else the quadratic step is taken.
        %
        info  = 2;
        bound = false;
        theta = 3*(LO.f - M.f)/(M.alpha - LO.alpha) + LO.Df + M.Df;
        s     = max(abs([theta,LO.Df,M.Df]));
        gamma = s*sqrt((theta/s)^2 - (LO.Df/s)*(M.Df/s));
        if M.alpha > LO.alpha
          gamma = -gamma;
        end
        p = (gamma - M.Df) + theta;
        q = ((gamma - M.Df) + gamma) + LO.Df;
        r = p/q;
        alphac = M.alpha + r*(LO.alpha - M.alpha);
        alphaq = M.alpha + (M.Df/(M.Df-LO.Df))*(LO.alpha - M.alpha);
        if abs(alphac-M.alpha) >  abs(alphaq-M.alpha)
          alphaf = alphac;
        else
          alphaf = alphaq;
        end
        bracketed = true;
      elseif abs(M.Df) < abs(LO.Df)
        %
        % Third case.
        % A lower function value, derivatives of the same sign, and the
        % magnitude of the derivative decreases.
        % The cubic step is only used if the cubic tends to infinity in the
        % direction of the step or if the minimum of the cubic is beyond M.alpha.
        % Otherwise the cubic step is defined to be either alphamin or alphamax.
        % The quadratic (secant) step is also computed and if the minimum is
        % bracketed then the the step closest to LO.alpha is taken,
        % else the step farthest away is taken.
        %
        info  = 3;
        bound = true;
        theta = 3*(LO.f - M.f)/(M.alpha - LO.alpha) + LO.Df + M.Df;
        s     = max(abs([theta,LO.Df,M.Df]));
        %
        % The case gamma = 0 only arises if the cubic does not tend
        % to infinity in the direction of the step.
        %
        gamma = s*sqrt(max(0,(theta/s)^2- (LO.Df/s)*(M.Df/s)));
        if M.alpha > LO.alpha
          gamma = -gamma;
        end
        p = (gamma - M.Df) + theta;
        q = (gamma + (LO.Df - M.Df)) + gamma;
        r = p/q
        if r < 0 && gamma ~= 0
          alphac = M.alpha + r*(LO.alpha - M.alpha);
        elseif M.alpha > LO.alpha
          alphac = alphamax;
        else
          alphac = alphamin;
        end
        alphaq = M.alpha + (M.Df/(M.Df-LO.Df))*(LO.alpha - M.alpha);
        if bracketed
          if abs(M.alpha-alphac) < abs(M.alpha-alphaq)
            alphaf = alphac;
          else
            alphaf = alphaq;
          end
        else
          if abs(M.alpha-alphac) > abs(M.alpha-alphaq)
            alphaf = alphac;
          else
            alphaf = alphaq;
          end
        end
      else
        %
        % Fourth case.
        % A lower function value, derivatives of the same sign, and the
        % magnitude of the derivative does not decrease.
        % If the minimum is not bracketed, the step is either alphamin or
        % alphamax, else the cubic step is taken.
        %
        info  = 4;
        bound = false;
        if bracketed
          theta = 3*(M.f - HI.f)/(HI.alpha - M.alpha) + HI.Df + M.Df;
          s     = max(abs([theta,HI.Df,M.Df]));
          gamma = s*sqrt((theta/s)**2 - (HI.Df/s)*(M.Df/s));
          if M.alpha > HI.alpha
            gamma = -gamma;
          end
          p = (gamma - M.Df) + theta;
          q = ((gamma - M.Df) + gamma) + HI.Df;
          r = p/q;
          alphac = M.alpha + r*(HI.alpha - M.alpha);
          alphaf = alphac;
        elseif M.alpha > LO.alpha
          alphaf = alphamax;
        else
          alphaf = alphamin;
        end
      end
      %
      % Update the interval of uncertainty.
      % This update does not depend on the new step or the case analysis above.
      %
      if M.f > LO.f
        HI = M;
      else
        if sgnd < 0
           HI = LO;
        end
        LO = M;
      end
      %
      % Compute the new step and safeguard it.
      %
      alphaf = min(alphamax,max(alphamin,alphaf);
      if bracketed && bound
        atmp = LO.alpha+0.66*(HI.alpha-LO.alpha)
        if HI.alpha > LO.alpha
          alphaf = min(atmp,alphaf);
        else
          alphaf = max(atmp,alphaf);
        end
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = LinesearchMoreThuente()
      self@LinesearchForwardBackward('MoreThuente');
      self.xtol         = 1e-2;
      self.xtrapf       = 4;
      self.max_fun_eval = 1000;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotDebug( self, alpha_guess )
      self.printInfo();
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function info = getInfo( self )
      info = self.info;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function printInfo( self )
      switch (self.info)
      case 0
        M.frintf('Improper input parameters.\n');
      case 1
        M.frintf('The sufficient decrease condition and the directional derivative condition hold.\n');
      case 2
        M.frintf('Relative width of the interval of uncertainty is at most xtol.\n');
      case 3
        M.frintf('Number of calls to fcn has reached max_fun_eval.\n');
      case 4
        M.frintf('The step is at the lower bound alpha_min.\n');
      case 5
        M.frintf('The step is at the upper bound alpha_max.\n');
      case 6
        M.frintf('Rounding errors prevent further progress.\n');
        M.frintf('There may not be a step which satisfies the\n');
        M.frintf('sufficient decrease and curvature conditions.\n');
        M.frintf('Tolerances may be too small..\n');
      end
    end
    %
    % ----------------------------------------------------------------------
    % Line Search Routine
    % ----------------------------------------------------------------------
    %
      subroutine cvsrch(n,x,f,g,s,stp,ftol,gtol,xtol,
     *                  stpmin,stpmax,maxfev,info,nfev,wa)
      integer n,maxfev,info,nfev
      double precision f,stp,ftol,gtol,xtol,stpmin,stpmax
      double precision x(n),g(n),s(n),wa(n)

      % The purpose of ******** is to find a step which satisfies
      % a sufficient decrease condition and a curvature condition.
      %
      % At each stage the subroutine updates an interval of uncertainty with
      % endpoints LO.alpha and HI.alpha. The interval of uncertainty is initially chosen
      % so that it contains a minimizer of the modified function
      %
      %     f(stp) - f(0) - ftol*stp*f'(0).
      %
      % If a step is obtained for which the modified function has a nonpositive
      % function value and nonnegative derivative, then the interval of
      % uncertainty is chosen so that it contains a minimizer of f(stp).
      %
      % The algorithm is designed to find a step which satisfies the sufficient
      % decrease condition
      %
      %     f(stp) <= f(0) + ftol*stp*f'(0),
      %
      % and the curvature condition
      %
      %     abs(f'(stp)) <=  gtol*abs(f'(0)).
      %
      % If ftol is less than gtol and if, for example, the function
      % is bounded below, then there is always a step which satisfies both
      % conditions. If no step can be found which satisfies both conditions,
      % then the algorithm usually stops when rounding errors prevent further
      % progress.
      % In this case stp only satisfies the sufficient decrease condition.
      %
      %	stp is a nonnegative variable.
      %     On input stp contains an initial estimate of a satisfactory step.
      %     On output stp contains the final estimate.
      %
      % ftol and gtol are nonnegative input variables.
      %      Termination occurs when the sufficient decrease condition and the
      %      directional derivative condition are satisfied.
      %
      %	xtol is a nonnegative input variable.
      %      Termination occurs  when the relative width of the interval
      %      of uncertainty is at most xtol.
      %
      %	stpmin and stpmax are nonnegative input variables which
      %	  specify lower and upper bounds for the step.
      %
      %	maxfev is a positive integer input variable.
      %        Termination occurs when the number of calls to fcn is at least
      %        maxfev by the end of an iteration.
      %
      %	info is an integer output variable set as follows:
      %
      %	  info = 0  Improper input parameters.
      %
      %	  info = 1  The sufficient decrease condition and the
      %             directional derivative condition hold.
      %
      %	  info = 2  Relative width of the interval of uncertainty
      %		        is at most xtol.
      %
      %	  info = 3  Number of calls to fcn has reached maxfev.
      %
      %	  info = 4  The step is at the lower bound stpmin.
      %
      %	  info = 5  The step is at the upper bound stpmax.
      %
      %	  info = 6  Rounding errors prevent further progress.
      %             There may not be a step which satisfies the
      %             sufficient decrease and curvature conditions.
      %             Tolerances may be too small.
      %
      % nfev is an integer output variable set to the number of calls to fcn.
      %
      xtrapf = 4;
      info   = 0;
      infoc  = 1;
      %
      % Check the input parameters for errors.
      %
      if stp <= 0 || ftol < 0 || gtol < 0 || xtol < 0 || stpmin < 0 || ...
         stpmax < stpmin || maxfev <= 0
        return;
      end
      %
      % Compute the initial gradient in the search direction
      % and check that s is a descent direction.
      %
      Df0 = dot(g,s);
      if Df0 <= 0
        return;
      end
      %
      % Initialize local variables.
      %
      bracketed = false;
      stage1    = true;
      nfev      = 0;
      f0        = f;
      c1Df0     = ftol*Df0;
      width     = stpmax - stpmin;
      width1    = width/0.5;
      %%      do 20 j = 1, n
      %%         wa(j) = x(j)
      %%  20    continue
      %
      % The variables LO.alpha, LO.f, LO.Df contain the values of the step,
      % function, and directional derivative at the best step.
      % The variables HI.alpha, HI.f, HI.Df contain the value of the step,
      % function, and derivative at the other endpoint of
      % the interval of uncertainty.
      % The variables stp, f, dg contain the values of the step,
      % function, and derivative at the current step.
      %
      LO.alpha = 0;
      LO.f     = f0;
      LO.Df    = Df0;
      HI       = LO;
      %
      % Start of iteration.
      %
      while true
        %
        % Set the minimum and maximum steps to correspond
        % to the present interval of uncertainty.
        %
        if bracketed
          stmin = min(LO.alpha,HI.alpha);
          stmax = max(LO.alpha,HI.alpha);
        else
          stmin = LO.alpha;
          stmax = stp + xtrapf*(stp - LO.alpha);
        end
        %
        % Force the step to be within the bounds stpmax and stpmin.
        %
        stp = max(stp,min(stp,stpmax));
        %
        % If an unusual termination is to occur then let
        % stp be the lowest point obtained so far.
        %
        if ((bracketed .and. (stp .le. stmin .or. stp .ge. stmax))
     *      .or. nfev .ge. maxfev-1 .or. infoc .eq. 0
     *      .or. (bracketed .and. stmax-stmin .le. xtol*stmax)) stp = LO.alpha
        %
        % Evaluate the function and gradient at stp
        % and compute the directional derivative.
        %
        [f,dg] = eval_FG(stp);
        ftest1 = f0 + stp*c1Df0;
        %
        % Test for convergence.
        %
        if (bracketed && (stp <= stmin || stp >= stmax) ) || infoc == 0
          info = 6;
        elseif stp == stpmax && f <= ftest1 && dg <= c1Df0
          info = 5;
        elseif stp == stpmin && (f > ftest1 || dg >= c1Df0)
          info = 4;
        elseif nfev >= maxfev
          info = 3;
        elseif bracketed && stmax-stmin <= xtol*stmax
          info = 2;
        elseif f <= ftest1 && abs(dg) <= gtol*(-Df0)
          info = 1;
        end
        %
        % Check for termination.
        %
        if info ~= 0
          return;
        end
        %
        % In the first stage we seek a step for which the modified
        % function has a nonpositive value and nonnegative derivative.
        %
        if stage1
          if f <= ftest1 && dg >= min(ftol,gtol)*Df0
            stage1 = false;
          end
        end
        %
        % A modified function is used to predict the step only if we have not
        % obtained a step for which the modified function has a nonpositive
        % function value and nonnegative derivative, and if a lower function
        % value has been obtained but the decrease is not sufficient.
        %
        if stage1 && f <= LO.f && f > ftest1
          %
          % Define the modified function and derivative values.
          %
          M.alpha = stp;
          M.f     = f - stp*c1Df0;
          M.Df    = dg - c1Df0;;

          LO.f  = LO.f - LO.alpha*c1Df0;
          LO.Df = LO.Df - c1Df0;

          HI.f  = HI.f - HI.alpha*c1Df0;
          HI.Df = HI.Df - c1Df0;
          %
          % Call cstep to update the interval of uncertainty
          % and to compute the new step.
          %
          [LO, HI, stp, bracketed, infoc] = self.safeguardedStep( LO, HI, M, bracketed );
          %
          % Reset the function and gradient values for f.
          %
          LO.f  = LO.f + LO.alpha*c1Df0
          LO.Df = LO.Df + c1Df0

          HI.f  = HI.f + HI.alpha*c1Df0
          HI.Df = HI.Df + c1Df0
        else
          %
          % Call cstep to update the interval of uncertainty
          % and to compute the new step.
          %
          M.alpha = stp;
          M.f     = f;
          M.Df    = dg;
          [LO, HI, stp, bracketed, infoc] = self.safeguardedStep( LO, HI, M, bracketed );
        end
        %
        % Force a sufficient decrease in the size of the interval of uncertainty.
        %
        if bracketed
          if abs(HI.alpha-LO.alpha) >= 0.66*width1
            stp = LO.alpha + 0.5*(HI.alpha - LO.alpha);
          end
          width1 = width;
          width  = abs(HI.alpha-LO.alpha);
        end
        %
        % End of iteration.
        %
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [alpha,ok] = search( self, alpha_guess )
      self.info  = 0;
      infoc      = 1;
      ok         = false;
      %
      % Check the input parameters for errors.
      %
      alpha = alpha_guess;
      if alpha <= 0; return; end

      [LO,HI,ierr] = self.ForwardBackward( alpha_guess );

      %
      % Compute the initial gradient in the search direction
      % and check that s is a descent direction.
      %
      fun  = @(alpha) self.fun1D.eval(alpha+LO.alpha);
      Dfun = @(alpha) self.fun1D.eval_D(alpha+LO.alpha);

      Df0 = Dfun(0);
      if Df0 >= 0; return; end
      %
      % Initialize local variables.
      %
      bracketed     = false;
      stage1     = true;
      n_fun_eval = 0;
      f0         = fun(0);
      c1Df0      = self.c1*Df0;
      width      = self.alpha_max - self.alpha_min;
      width1     = 2*width;
      %
      % The variables LO.alpha, LO.f, DLO.f contain the values of the step,
      % function, and directional derivative at the best step.
      % The variables HI.alpha, HI.f, DHI.f contain the value of the step,
      % function, and derivative at the other enM.Dfoint of
      % the interval of uncertainty.
      % The variables alpha, f, Df contain the values of the step,
      % function, and derivative at the current step.
      %
      LO.alpha = 0.0;
      LO.f  = f0;
      DLO.f = Df0;
      HI.alpha = 0.0;
      HI.f  = f0;
      DHI.f = Df0;
      %
      % Start of iteration.
      %
      while true
        %
        % Set the minimum and maximum steps to correspond
        % to the present interval of uncertainty.
        %
        if bracketed
          stmin = min(LO.alpha,HI.alpha);
          stmax = max(LO.alpha,HI.alpha);
        else
          stmin = LO.alpha;
          stmax = alpha + self.xtrapf*(alpha - LO.alpha);
        end
        %
        % Force the step to be within the bounds alpha_max and alpha_min.
        %
        alpha = min(max(alpha,self.alpha_min),self.alpha_max);
        %
        % If an unusual termination is to occur then let
        % alpha be the lowest point obtained so far.
        %
        if ( bracketed && ( alpha <= stmin || ...
                         alpha >= stmax || ...
                         stmax-stmin <= self.xtol*stmax ) ) ...
           || n_fun_eval >= self.max_fun_eval || infoc == 0
          alpha = LO.alpha;
        end
        %
        % Evaluate the function and gradient at alpha
        % and compute the directional derivative.
        %
        f  = fun(alpha);
        Df = Dfun(alpha);
        n_fun_eval = n_fun_eval + 1;
        ftest1 = f0 + alpha*c1Df0;
        %
        % Test for convergence.
        %
        if (bracketed && (alpha <= stmin || alpha >= stmax)) || infoc == 0
          self.info = 6;
        elseif alpha == self.alpha_max && f <= ftest1 && Df <= c1Df0
          self.info = 5;
        elseif alpha == self.alpha_min && ( f > ftest1 || Df >= c1Df0 )
          self.info = 4;
        elseif n_fun_eval >= self.max_fun_eval
          self.info = 3;
        elseif bracketed && stmax-stmin <= self.xtol*stmax
          self.info = 2;
        elseif f <= ftest1 && abs(Df) <= self.c2*(-Df0)
          self.info = 1;
        end
        %
        % Check for termination.
        %
        %if self.info ~= 0; ok = self.info == 1; return; end
        if self.info ~= 0; ok = true; alpha = alpha + LO.alpha; return; end
        %
        % In the first stage we seek a step for which the modified
        % function has a nonpositive value and nonnegative derivative.
        %
        if stage1
          if f <= ftest1 && Df >= min(self.c1,self.c2)*Df0
            stage1 = false;
          end
        end
        %
        % A modified function is used to predict the step only if we have not
        % obtained a step for which the modified function has a nonpositive
        % function value and nonnegative derivative, and if a lower function
        % value has been obtained but the decrease is not sufficient.
        %
        if stage1 && f <= LO.f && f > ftest1
          %
          % Define the modified function and derivative values.
          %
          fm   = f - alpha*c1Df0;
          LO.fm  = LO.f - LO.alpha*c1Df0;
          HI.fm  = HI.f - HI.alpha*c1Df0;
          Dfm  = Df - c1Df0;
          DLO.fm = DLO.f - c1Df0;
          DHI.fm = DHI.f - c1Df0;
          %
          % Call cstep to update the interval of uncertainty
          % and to compute the new step.
          %
          [LO.alpha,LO.fm,DLO.fm,HI.alpha,HI.fm,DHI.fm,alpha,fm,Dfm,bracketed,infoc] ...
            = self.safeguardedStep(LO.alpha,LO.fm,DLO.fm,HI.alpha,HI.fm,DHI.fm,alpha,fm,Dfm,bracketed,stmin,stmax);
          %
          % Reset the function and gradient values for f.
          %
          LO.f  = LO.fm + LO.alpha*c1Df0;
          HI.f  = HI.fm + HI.alpha*c1Df0;
          DLO.f = DLO.fm + c1Df0;
          DHI.f = DHI.fm + c1Df0;
        else
          %
          % Call cstep to update the interval of uncertainty
          % and to compute the new step.
          %
          [LO.alpha,LO.f,DLO.f,HI.alpha,HI.f,DHI.f,alpha,f,Df,bracketed,infoc] ...
            = self.safeguardedStep(LO.alpha,LO.f,DLO.f,HI.alpha,HI.f,DHI.f,alpha,f,Df,bracketed,stmin,stmax);
        end
        %
        % Force a sufficient decrease in the size of the
        % interval of uncertainty.
        %
        if bracketed
          if abs(HI.alpha-LO.alpha) >= 0.66*width1
            alpha = LO.alpha + 0.5*(HI.alpha - LO.alpha);
          end
          width1 = width;
          width  = abs(HI.alpha-LO.alpha);
        end
        % End of iteration.
      end
    end
  end
end
