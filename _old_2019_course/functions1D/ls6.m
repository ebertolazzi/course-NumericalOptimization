classdef ls6 < Function1D
  properties (SetAccess = private, Hidden = true)
  end
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = ls6()
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function r = getRange( self )
      r = [0,1];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, alpha )
      tmp1 = (alpha-0.6)./ 0.4;
      tmp2 = (alpha-0.2)./ 0.04;
      y = 2 - 0.8*exp(-tmp1.^2 ) - exp( -tmp2.^2 ) ;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, alpha )
      tmp1  = (alpha-0.6)./ 0.4;
      tmp2  = (alpha-0.2)./ 0.04;
      Dtmp1 = 1/0.4;
      Dtmp2 = 1/0.04;
      Dy   = 0.8*(2*tmp1*Dtmp1).* exp(-tmp1.^2 ) + (2*tmp2*Dtmp2) .* exp( -tmp2.^2 );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, alpha )
      tmp1  = (alpha-0.6)./ 0.4;
      tmp2  = (alpha-0.2)./ 0.04;
      Dtmp1 = 1/0.4;
      Dtmp2 = 1/0.04;
      DDy   = 0.8*(2*Dtmp1^2).* exp(-tmp1.^2 ) + (2*Dtmp2^2) .* exp( -tmp2.^2 ) + ...
            - 0.8*(2*tmp1*Dtmp1).^2.* exp(-tmp1.^2 ) - (2*tmp2*Dtmp2).^2 .* exp( -tmp2.^2 );
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
