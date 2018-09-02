classdef LinesearchArmijo < LinesearchForwardBackward
  % Armijo linesearch

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = LinesearchArmijo()
      self@LinesearchForwardBackward('Armijo');
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [alpha_ott,ok] = search( self, alpha_guess )
      % find step that satisfy Armijo condition
      % of ok = false search failed
      [alpha0,alpha1,~,~,ierr] = self.ForwardBackward( alpha_guess );
      ok = true;
      switch ierr
      case 0 % only alpha0 satisfy Armijo
        alpha_ott = alpha0;
      case {1,2} % only alpha0 and alpha1 satisfy Armijo, take average
        alpha_ott = (alpha0+alpha1)/2;
      otherwise
        alpha_ott = 0;
        ok        = false;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
