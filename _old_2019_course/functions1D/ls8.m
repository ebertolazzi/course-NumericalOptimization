classdef ls8 < Function1D
  properties (SetAccess = private, Hidden = true)
  end
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = ls8()
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function r = getRange( self )
      r = [0,0.4];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, alpha )
      idx    = find( alpha > 0);
      a      = alpha(idx);
      y(idx) = -a + 1e3*a.^3.*sin(1./a);
      idx    = find( alpha <= 0);
      y(idx) = 0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, alpha )
      idx     = find( alpha > 0);
      a       = alpha(idx);
      ra      = 1./a;
      Dy(idx) = -1 + 1e3*a.*(3*a*sin(ra)-cos(ra));
      idx     = find( alpha <= 0);
      Dy(idx) = 0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, alpha )
      idx     = find( alpha > 0);
      ra      = 1./alpha(idx);
      Dy(idx) = 1e3*(6*alpha(idx)-ra).*sin(ra) - 4e3*cos(ra);
      idx     = find( alpha <= 0);
      Dy(idx) = 0;
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
