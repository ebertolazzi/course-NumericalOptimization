classdef testFun3 < Function1D
  properties (SetAccess = private, Hidden = true)
    epsi
    xsing
  end
  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = testFun3()
      self.epsi  = 0.1 ;
      self.xsing = 5;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, xx )
      y      = Inf*ones(size(xx));
      idx    = find( xx < self.xsing );
      x      = xx(idx);
      y(idx) = exp(-3*x)+x.^2.*sin(3*x)+self.epsi./(self.xsing-x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, xx )
      Dy      = Inf*ones(size(xx));
      idx     = find( xx < self.xsing );
      x       = xx(idx);
      Dy(idx) = -3*exp(-3*x)+2*x*sin(3*x)+3*x^2*cos(3*x)+self.epsi./(self.xsing-x).^2;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, xx )
      DDy      = Inf*ones(size(xx));
      idx      = find( xx < self.xsing );
      x        = xx(idx);
      DDy(idx) = 9*exp(-3*x)+(2-9*x^2)*sin(3*x)+12*x*cos(3*x)+2*self.epsi./(self.xsing-x).^3;
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
