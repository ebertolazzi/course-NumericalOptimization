classdef dummy < handle
  properties (SetAccess = private, Hidden = true)
    epsi;
  end
  
  methods (Access = private)
    function res = sum2( self, x, y )
      res = x^2+y^2;
    end
  end
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = dummy()
      self.epsi = 1e-8;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,varargin)
      if nargin > 2
        disp('use mode 1');
        x = varargin{1};
        y = varargin{2};
        f = self.sum2(x,y);
      else
        disp('use mode 2');
        xy = varargin{1};
        f = self.sum2(xy(1),xy(1));
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
