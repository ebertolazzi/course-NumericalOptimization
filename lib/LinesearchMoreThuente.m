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
  %  uncertainty with endpoints stx and sty.
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
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [stx,fx,Dfx,sty,fy,Dfy,stp,fp,Dfp,brackt,info] ...
      = safeguardedStep(self,stx,fx,Dfx,sty,fy,Dfy,stp,fp,Dfp,brackt,stpmin,stpmax)
      %%
      % The purpose of cstep is to compute a safeguarded step for
      % a linesearch and to update an interval of uncertainty for
      % a minimizer of the function.
      %
      % The parameter stx contains the step with the least function value.
      % The parameter stp contains the current step.
      % It is assumed that the derivative at stx is negative in the
      % direction of the step. If brackt is set true then a minimizer has been
      % bracketed in an interval of uncertainty with endpoints stx and sty.
      %
      % stx, fx, and Dfx are variables which specify the step, the function,
      % and the derivative at the best step obtained so far.
      % The derivative must be negative in the direction of the step, that is,
      % Dfx and stp-stx must have opposite signs.
      % On output these parameters are updated appropriately.
      %
      % sty, fy, and Dfy are variables which specify the step, the function,
      % and the derivative at the other endpoint of the interval of uncertainty.
      % On output these parameters are updated appropriately.
      %
      % stp, fp, and Dfp are variables which specify the step, the function,
      % and the derivative at the current step.
      % If brackt is set true then on input stp must be between stx and sty.
      % On output stp is set to the new step.
      %
      % brackt is a logical variable which specifies if a minimizer
      % has been bracketed. If the minimizer has not been bracketed
      % then on input brackt must be set false. If the minimizer
      % is bracketed then on output brackt is set true.
      %
      % stpmin and stpmax are input variables which specify lower
      % and upper bounds for the step.
      %
      % info is an integer output variable set as follows:
      % If info = 1,2,3,4,5, then the step has been computed
      % according to one of the five cases below. Otherwise
      % info = 0, and this indicates improper input parameters.
      %
      info = 0;
      %
      %  Check the input parameters for errors.
      %
      if (brackt && (stp <= min(stx,sty) || stp >= max(stx,sty))) ...
        || Dfx*(stp-stx) >= 0.0 || stpmax < stpmin
        return;
      end
      %
      % Determine if the derivatives have opposite sign.
      %
      sgnd = Dfp*sign(Dfx);
      %
      % First case. A higher function value. The minimum is bracketed.
      % If the cubic step is closer to stx than the quadratic step,
      % the cubic step is taken, else the average of the cubic and
      % quadratic steps is taken.
      %
      if fp > fx
        info  = 1;
        bound = 1;
        theta = 3*(fx - fp)/(stp - stx) + Dfx + Dfp;
        s     = norm([theta,Dfx,Dfp],inf);
        gamma = s*sqrt((theta/s)^2 - (Dfx/s)*(Dfp/s));
        if (stp < stx); gamma = -gamma; end
        p = (gamma - Dfx) + theta;
        q = ((gamma - Dfx) + gamma) + Dfp;
        r = p/q;
        stpc = stx + r*(stp - stx);
        stpq = stx + ((Dfx/((fx-fp)/(stp-stx)+Dfx))/2)*(stp - stx);
        if abs(stpc-stx) < abs(stpq-stx)
          stpf = stpc;
        else
          stpf = stpc + (stpq - stpc)/2;
        end
        brackt = true;
        %
        % Second case. A lower function value and derivatives of opposite sign.
        % The minimum is bracketed. If the cubic step is closer to stx than
        % the quadratic (secant) step, the cubic step is taken,
        % else the quadratic step is taken.
        %
      elseif sgnd < 0
        info  = 2;
        bound = false;
        theta = 3*(fx - fp)/(stp - stx) + Dfx + Dfp;
        s     = norm([theta,Dfx,Dfp],inf);
        gamma = s*sqrt((theta/s)^2 - (Dfx/s)*(Dfp/s));
        if stp > stx; gamma = -gamma; end
        p    = (gamma - Dfp) + theta;
        q    = ((gamma - Dfp) + gamma) + Dfx;
        r    = p/q;
        stpc = stp + r*(stx - stp);
        stpq = stp + (Dfp/(Dfp-Dfx))*(stx - stp);
        if abs(stpc-stp) > abs(stpq-stp)
          stpf = stpc;
        else
          stpf = stpq;
        end
        brackt = true;
        %
        % Third case. A lower function value, derivatives of the same sign,
        % and the magnitude of the derivative decreases.
        % The cubic step is only used if the cubic tends to infinity in the
        % direction of the step or if the minimum of the cubic is beyond stp.
        % Otherwise the cubic step is defined to be either stpmin or stpmax.
        % The quadratic (secant) step is also computed and if the minimum is
        % bracketed then the the step closest to stx is taken, else the step
        % farthest away is taken.
        %
      elseif abs(Dfp) < abs(Dfx)
        info  = 3;
        bound = 1;
        theta = 3*(fx - fp)/(stp - stx) + Dfx + Dfp;
        s     = norm([theta,Dfx,Dfp],inf);
        %
        % The case gamma = 0 only arises if the cubic does not tend to
        % infinity in the direction of the step.
        %
        gamma = s*sqrt(max(0.,(theta/s)^2 - (Dfx/s)*(Dfp/s)));
        if stp > stx; gamma = -gamma; end
        p = (gamma - Dfp) + theta;
        q = (gamma + (Dfx - Dfp)) + gamma;
        r = p/q;
        if r < 0.0 && gamma ~= 0
          stpc = stp + r*(stx - stp);
        elseif stp > stx
          stpc = stpmax;
        else
          stpc = stpmin;
        end
        stpq = stp + (Dfp/(Dfp-Dfx))*(stx - stp);
        if brackt
          if abs(stp-stpc) < abs(stp-stpq)
            stpf = stpc;
          else
            stpf = stpq;
          end
        else
          if abs(stp-stpc) > abs(stp-stpq)
            stpf = stpc;
          else
            stpf = stpq;
          end
        end
        %
        % Fourth case. A lower function value, derivatives of the same sign,
        % and the magnitude of the derivative does not decrease.
        % If the minimum is not bracketed, the step is either stpmin or stpmax,
        % else the cubic step is taken.
        %
      else
        info  = 4;
        bound = false;
        if brackt
          theta = 3*(fp - fy)/(sty - stp) + Dfy + Dfp;
          s     = norm([theta,Dfy,Dfp],inf);
          gamma = s*sqrt((theta/s)^2 - (Dfy/s)*(Dfp/s));
          if stp > sty; gamma = -gamma; end
          p    = (gamma - Dfp) + theta;
          q    = ((gamma - Dfp) + gamma) + Dfy;
          r    = p/q;
          stpc = stp + r*(sty - stp);
          stpf = stpc;
        elseif stp > stx
          stpf = stpmax;
        else
          stpf = stpmin;
        end
      end
      %
      % Update the interval of uncertainty.
      % This update does not depend on the new step or the case analysis above.
      %
      if fp > fx
        sty = stp;
        fy  = fp;
        Dfy = Dfp;
      else
        if sgnd < 0.0
          sty = stx;
          fy  = fx;
          Dfy = Dfx;
        end
        stx = stp;
        fx  = fp;
        Dfx = Dfp;
      end
      %
      % Compute the new step and safeguard it.
      %
      stpf = min(stpmax,stpf);
      stpf = max(stpmin,stpf);
      stp  = stpf;
      if brackt && bound
        if sty > stx
          stp = min(stx+0.66*(sty-stx),stp);
        else
          stp = max(stx+0.66*(sty-stx),stp);
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
        fprintf('Improper input parameters.\n');
      case 1
        fprintf('The sufficient decrease condition and the directional derivative condition hold.\n');
      case 2
        fprintf('Relative width of the interval of uncertainty is at most xtol.\n');
      case 3
        fprintf('Number of calls to fcn has reached max_fun_eval.\n');
      case 4
        fprintf('The step is at the lower bound alpha_min.\n');
      case 5
        fprintf('The step is at the upper bound alpha_max.\n');
      case 6
        fprintf('Rounding errors prevent further progress.\n');
        fprintf('There may not be a step which satisfies the\n');
        fprintf('sufficient decrease and curvature conditions.\n');
        fprintf('Tolerances may be too small..\n');
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
      brackt     = false;
      stage1     = true;
      n_fun_eval = 0;
      f0         = fun(0);
      c1Df0      = self.c1*Df0;
      width      = self.alpha_max - self.alpha_min;
      width1     = 2*width;
      %
      % The variables stx, fx, Dfx contain the values of the step,
      % function, and directional derivative at the best step.
      % The variables sty, fy, Dfy contain the value of the step,
      % function, and derivative at the other endpoint of
      % the interval of uncertainty.
      % The variables alpha, f, Df contain the values of the step,
      % function, and derivative at the current step.
      %
      stx = 0.0;
      fx  = f0;
      Dfx = Df0;
      sty = 0.0;
      fy  = f0;
      Dfy = Df0;
      %
      % Start of iteration.
      %
      while true
        %
        % Set the minimum and maximum steps to correspond
        % to the present interval of uncertainty.
        %
        if brackt
          stmin = min(stx,sty);
          stmax = max(stx,sty);
        else
          stmin = stx;
          stmax = alpha + self.xtrapf*(alpha - stx);
        end
        %
        % Force the step to be within the bounds alpha_max and alpha_min.
        %
        alpha = min(max(alpha,self.alpha_min),self.alpha_max);
        %
        % If an unusual termination is to occur then let
        % alpha be the lowest point obtained so far.
        %
        if ( brackt && ( alpha <= stmin || ...
                         alpha >= stmax || ...
                         stmax-stmin <= self.xtol*stmax ) ) ...
           || n_fun_eval >= self.max_fun_eval || infoc == 0
          alpha = stx;
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
        if (brackt && (alpha <= stmin || alpha >= stmax)) || infoc == 0
          self.info = 6;
        elseif alpha == self.alpha_max && f <= ftest1 && Df <= c1Df0
          self.info = 5;
        elseif alpha == self.alpha_min && ( f > ftest1 || Df >= c1Df0 )
          self.info = 4;
        elseif n_fun_eval >= self.max_fun_eval
          self.info = 3;
        elseif brackt && stmax-stmin <= self.xtol*stmax
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
        if stage1 && f <= fx && f > ftest1
          %
          % Define the modified function and derivative values.
          %
          fm   = f - alpha*c1Df0;
          fxm  = fx - stx*c1Df0;
          fym  = fy - sty*c1Df0;
          Dfm  = Df - c1Df0;
          Dfxm = Dfx - c1Df0;
          Dfym = Dfy - c1Df0;
          %
          % Call cstep to update the interval of uncertainty
          % and to compute the new step.
          %
          [stx,fxm,Dfxm,sty,fym,Dfym,alpha,fm,Dfm,brackt,infoc] ...
            = self.safeguardedStep(stx,fxm,Dfxm,sty,fym,Dfym,alpha,fm,Dfm,brackt,stmin,stmax);
          %
          % Reset the function and gradient values for f.
          %
          fx  = fxm + stx*c1Df0;
          fy  = fym + sty*c1Df0;
          Dfx = Dfxm + c1Df0;
          Dfy = Dfym + c1Df0;
        else
          %
          % Call cstep to update the interval of uncertainty
          % and to compute the new step.
          %
          [stx,fx,Dfx,sty,fy,Dfy,alpha,f,Df,brackt,infoc] ...
            = self.safeguardedStep(stx,fx,Dfx,sty,fy,Dfy,alpha,f,Df,brackt,stmin,stmax);
        end
        %
        % Force a sufficient decrease in the size of the
        % interval of uncertainty.
        %
        if brackt
          if abs(sty-stx) >= 0.66*width1
            alpha = stx + 0.5*(sty - stx);
          end
          width1 = width;
          width  = abs(sty-stx);
        end
        % End of iteration.
      end
    end
  end
end
