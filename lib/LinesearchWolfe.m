classdef LinesearchWolfe < LinesearchForwardBackward
  % Wolfe linesearch

  properties (SetAccess = private, Hidden = true)
    c2      % Wolfe constant to accept the step (0,1/2]
            % f'(alpha) >= c2*f'(0)
    alpha_epsi
    dumpMin
    dumpMax
    strongWolfe
  end

  methods

    function self = LinesearchWolfe()
      self@LinesearchForwardBackward();
      self.c2          = max(self.c1,0.1) ;
      self.alpha_epsi  = 1e-10 ;
      self.dumpMin     = 0.95 ;
      self.dumpMax     = 0.05 ;
      self.strongWolfe = false ;
    end

    function setC1C2( self, c1, c2 )
      if c1 > c2
        error('LinesearchWolfe, constant c1 = %g must be <= c2 = %g\n',c1,c2);
      end
      self.setC1( c1 ) ;
      % set the c1 coefficients for lineasearch
      if c2 >= 1
        error('LinesearchWolfe, constant c2 = %g must be < 1\n',c2);
      end
      self.c2 = c2 ;
    end

    function setDump( self, d1, d2 )
      if d1 > 0.4
        error('LinesearchWolfe, dump d1 = %g must be <= 0.4%g\n',d1);
      end
      if d2 < 0.6
        error('LinesearchWolfe, dump d2 = %g must be >= 0.6%g\n',d2);
      end
      self.dumpMin = d1 ;
      self.dumpMax = d2 ;
    end

    function strongWolfe_on( self )
      self.strongWolfe = true ;
    end

    function strongWolfe_off( self )
      self.strongWolfe = false ;
    end
    
    %                         _           _   _
    %    __ _ _   _  __ _  __| |_ __ __ _| |_(_) ___ 
    %   / _` | | | |/ _` |/ _` | '__/ _` | __| |/ __|
    %  | (_| | |_| | (_| | (_| | | | (_| | |_| | (__ 
    %   \__, |\__,_|\__,_|\__,_|_|  \__,_|\__|_|\___|
    %      |_|                                       
    %      
    function alpha = quadratic( ~, f0, Df0, fp, p )
      alpha = Df0 * p^2 / ( 2*(f0+Df0*p-fp) ) ;
    end

    function [alpha_ott,ok] = search( self, alpha_guess )
      % find step that satisfy Armijo condition
      % of ok = false search failed
      [lLo,lHi,ierr] = self.ForwardBackward( alpha_guess ) ;
      
      % check if search is failed
      if ierr < 0
        ok = false ;
        return ;
      end

      f0    = self.fun1D.eval(0) ;
      Df0   = self.fun1D.eval_D(0) ;
      c1Df0 = self.c1 * Df0 ;
      c2Df0 = self.c2 * Df0 ;
      
      % check if lLo satisfy wolfe condition
      DfLo = self.fun1D.eval_D(lLo) ;
      if DfLo >= c2Df0 && ( ~self.strongWolfe || DfLo <= -c2Df0 )
        alpha_ott = lLo ;
        ok        = true ;
        return ;
      end

      % check if interval is left or right
      if DfLo > 0
        lHi = 0 ;
      end
      
      %   ____  _____ _____ ___ _   _ _____ 
      %  |  _ \| ____|  ___|_ _| \ | | ____|
      %  | |_) |  _| | |_   | ||  \| |  _|  
      %  |  _ <| |___|  _|  | || |\  | |___ 
      %  |_| \_\_____|_|   |___|_| \_|_____|
      %
      fLo = self.fun1D.eval(lLo) ;
      fHi = self.fun1D.eval(lHi) ;
 
      Delta = lHi - lLo ;
 
      while abs(Delta) > self.alpha_epsi

        deltaLambda = self.quadratic( fLo, DfLo, fHi, Delta ) ;
        if Delta > 0
          deltaLambda = min( Delta*self.dumpMax, max( Delta*self.dumpMin, deltaLambda )) ;
        else
          deltaLambda = max( Delta*self.dumpMax, min( Delta*self.dumpMin, deltaLambda )) ;
        end
       
        l  = lLo + deltaLambda ;
        fl = self.fun1D.eval(l) ;

        if fl > f0 + l*c1Df0 || fl > fLo
          lHi = l ;
          fHi = fl ;
        else
          Dfl = self.fun1D.eval_D(l) ;

          if Dfl >= c2Df0 && ( ~self.strongWolfe || Dfl <= -c2Df0 )
            alpha_ott = l ;
            ok        = true ;
            return ; % found Wolfe point
          end
          lLo  = l ;
          fLo  = fl ;
          DfLo = Dfl ;
        end
        Delta = lHi - lLo ;
      end
      error( 'WolfeLineSearch (refine): failed refine lHi=%g lLo=%g DfLo=%g\n', lHi, lLo, DfLo ) ;
    end
  end
end
