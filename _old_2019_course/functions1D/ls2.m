classdef ls2 < Function1D
  properties (SetAccess = private, Hidden = true)
    beta
  end
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = ls2()
      self.beta = 0.004;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function r = getRange( self )
      r = [0,2];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, alpha )
      y = (alpha+self.beta).^5 - 2*(alpha+self.beta).^4 ;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, alpha )
      Dy = 5*(alpha+self.beta).^4 - 8*(alpha+self.beta).^3 ;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, alpha )
      Dy = 20*(alpha+self.beta).^3 - 24*(alpha+self.beta).^2 ;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [y,Dy] = eval_FG( self, x )
      y  = self.eval(x);
      Dy = self.eval_D(x);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [y,Dy,DDy] = eval_FGH( self, x )
      y   = self.eval(x);
      Dy  = self.eval_D(x);
      DDy = self.eval_DD(x);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
