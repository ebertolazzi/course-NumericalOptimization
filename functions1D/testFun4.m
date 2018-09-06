classdef testFun4 < Function1D
  properties (SetAccess = private, Hidden = true)
    nnn
    slope
    epsi
    xsing
  end
  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = testFun4()
      self.nnn   = 50;
      self.slope = 0.2;
      self.epsi  = 0.1;
      self.xsing = 5;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, xx )
      y      = Inf*ones(size(xx));
      idx    = find( xx >= 0 & xx < self.xsing );
      x      = xx(idx);
      y(idx) = self.slope*x + ...
               1./(1+self.nnn*sqrt(x)) + ...
               self.epsi./(self.xsing-x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, xx )
      Dy      = -Inf*ones(size(xx));
      idx     = find( xx >= 0 & xx < self.xsing  );
      x       = xx(idx);
      Dy(idx) = self.slope - ...
                self.nnn./(2*sqrt(x).*(1+self.nnn*sqrt(x)).^2) + ...
                self.epsi./(self.xsing-x).^2;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, xx )
      DDy      = Inf*ones(size(xx));
      idx      = find( xx >= 0 & xx < self.xsing );
      x        = xx(idx);
      DDy(idx) = (self.nnn/4).*(3*self.nnn*sqrt(x)+1)/(x.^(3/2).*(1+self.nnn*sqrt(x)).^3) + ...
                 2*self.epsi./(self.xsing-x).^3;
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
