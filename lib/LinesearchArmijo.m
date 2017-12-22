classdef LinesearchArmijo < LinesearchForwardBackward
  % Armijo linesearch

  methods

    function self = LinesearchArmijo()
      self@LinesearchForwardBackward();
    end

    function [alpha_ott,ok] = search( self, alpha_guess )
      % find step that satisfy Armijo condition
      % of ok = false search failed
      [alpha0,~,ierr] = self.ForwardBackward( alpha_guess );
      alpha_ott = alpha0 ;
      ok        = ierr >= 0 ;
    end

  end
end
