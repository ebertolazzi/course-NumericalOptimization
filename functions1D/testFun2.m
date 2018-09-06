classdef testFun2 < Function1D
  properties (SetAccess = private, Hidden = true)
    epsi
    epsi2
    xsing
  end
  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = testFun2()
      self.epsi  = 0.1;
      self.epsi2 = 10;
      self.xsing = 5;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, xx )
      y      = Inf*ones(size(xx));
      idx    = find( xx < self.xsing );
      x      = xx(idx);
      y(idx) = -self.epsi2./(1+(x-1).^2)-x.^2+self.epsi./(self.xsing-x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, xx )
      Dy      = Inf*ones(size(xx));
      idx     = find( xx < self.xsing );
      x       = xx(idx);
      Dy(idx) = -2*self.epsi2*(1-x)./(1+(x-1).^2).^2-2*x+self.epsi./(self.xsing-x).^2;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, xx )
      DDy      = Inf*ones(size(xx));
      idx      = find( xx < self.xsing );
      x        = xx(idx);
      DDy(idx) = -self.epsi2*(6*x.^2-12*x+4)./(1+(x-1).^2).^3-2+2*self.epsi./(self.xsing-x).^3;
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
