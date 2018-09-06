classdef testFun1 < Function1D
  properties (SetAccess = private, Hidden = true)
    epsi
    xsing
  end
  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = testFun1()
      self.epsi  = 0.1 ;
      self.xsing = 5;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, xx )
      y      = Inf*ones(size(xx));
      idx    = find( xx < self.xsing );
      x      = xx(idx);
      y(idx) = -x+self.epsi./(self.xsing-x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, xx )
      Dy      = Inf*ones(size(xx));
      idx     = find( xx < self.xsing );
      x       = xx(idx);
      Dy(idx) = -1+self.epsi./(self.xsing-x).^2;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, xx )
      DDy      = Inf*ones(size(xx));
      idx      = find( xx < self.xsing );
      x        = xx(idx);
      DDy(idx) = 2*self.epsi./(self.xsing-x).^3;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [y,Dy] = eval_FG( self, x )
      y  = self.eval(x);
      Dy = self.eval_D(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [y,Dy,DDy] = eval_FGH( self, x )
      y   = self.eval(x);
      Dy  = self.eval_D(x);
      DDy = self.eval_DD(x);
    end
  end
end
