classdef LinesearchWolfe < LinesearchForwardBackward
  % Wolfe linesearch

  properties (SetAccess = private, Hidden = true)
    strongWolfe
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = LinesearchWolfe()
      self@LinesearchForwardBackward('Wolfe');
      self.strongWolfe = false;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function strongWolfe_on( self )
      self.strongWolfe = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function strongWolfe_off( self )
      self.strongWolfe = false;
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
    function [alpha_ott,ok] = search( self, alpha_guess )
      % find step that satisfy Armijo condition
      % of ok = false search failed
      [LO,HI,ierr] = self.ForwardBackward( alpha_guess );
      switch ierr
      case 0 % only alpha0 satisfy Armijo, check if minimum is on [0,aLO]
        if LO.Df > 0.0
          [alpha_ott,ok] = self.Zoom( 0.0, self.f0, LO.alpha, LO.f, self.strongWolfe );
        else
          [alpha_ott,ok] = self.Zoom( LO.alpha, LO.f, HI.alpha, HI.f, self.strongWolfe );
        end
      %case 1 %
      %  [alpha_ott,ok] = self.Zoom( aLO, fLO, aHI, fHI, self.strongWolfe );
      %case 2 % LO and HI exchanged
      %  [alpha_ott,ok] = self.Zoom( aHI, fHI, aLO, fLO, self.strongWolfe );
      otherwise
        ok        = false;
        alpha_ott = 0;
      end
    end
  end
end
