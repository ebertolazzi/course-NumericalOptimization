classdef ls3 < Function1D
  properties (SetAccess = private, Hidden = true)
    beta
  end
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = ls3()
      self.beta = 100;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function r = getRange( self )
      r = [0,3];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, alpha )
      idx    = find( alpha < 0 );
      y(idx) = self.beta*alpha(idx).*(self.beta*alpha(idx)-1) ;
      idx    = find( alpha >= 0 & alpha < 1 );
      y(idx) = -log( 1+self.beta*alpha(idx));
      idx    = find( alpha >= 1 );
      t1     = self.beta/(1+self.beta);
      y(idx) = -log( 1+self.beta) + t1 *(alpha(idx)-1).*(t1 *(alpha(idx)-1)-1);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, alpha )
      idx     = find( alpha < 0 );
      Dy(idx) = self.beta*(2*self.beta*alpha(idx)-1) ;
      idx     = find( alpha >= 0 & alpha < 1 );
      Dy(idx) = -self.beta./(1+self.beta*alpha(idx));
      idx     = find( alpha >= 1 );
      t1      = self.beta/(1+self.beta);
      Dy(idx) = t1*(2*t1 *(alpha(idx)-1)-1);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, alpha )
      idx      = find( alpha < 0 );
      DDy(idx) = 2*self.beta^2;
      idx      = find( alpha >= 0 & alpha < 1 );
      DDy(idx) = self.beta^2./(1+self.beta*alpha(idx)).^2;
      idx      = find( alpha >= 1 );
      t1       = self.beta/(1+self.beta);
      DDy(idx) = 2*t1^2;
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
