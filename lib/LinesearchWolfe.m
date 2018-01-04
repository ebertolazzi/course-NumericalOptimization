classdef LinesearchWolfe < LinesearchForwardBackward
  % Wolfe linesearch

  properties (SetAccess = private, Hidden = true)
    strongWolfe
  end

  methods

    function self = LinesearchWolfe()
      self@LinesearchForwardBackward('Wolfe');
      self.strongWolfe = false ;
    end

    function strongWolfe_on( self )
      self.strongWolfe = true ;
    end

    function strongWolfe_off( self )
      self.strongWolfe = false ;
    end

    function [alpha_ott,ok] = search( self, alpha_guess )
      % find step that satisfy Armijo condition
      % of ok = false search failed
      [aLO,aHI,fLO,fHI,ierr] = self.ForwardBackward( alpha_guess ) ;
      switch ierr
      case 0 % only alpha0 satisfy Armijo, check if minimum is on [0,aLO]
        if self.fun1D.eval_D(aLO) > 0
          [alpha_ott,ok] = self.Zoom( 0, self.f0, aLO, fLO, self.strongWolfe ) ;
        else
          [alpha_ott,ok] = self.Zoom( aLO, fLO, aHI, fHI, self.strongWolfe ) ;
        end
      case 1 %
        [alpha_ott,ok] = self.Zoom( aLO, fLO, aHI, fHI, self.strongWolfe ) ;
      case 2 % LO and HI exchanged
        [alpha_ott,ok] = self.Zoom( aHI, fHI, aLO, fLO, self.strongWolfe ) ;
      otherwise
        ok        = false ;
        alpha_ott = 0 ;
      end
    end
  end
end
