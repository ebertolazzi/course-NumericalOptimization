classdef ls5 < Function1D
  properties (SetAccess = private, Hidden = true)
    beta1
    beta2
  end
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = ls5()
      self.beta1 = 0.001;
      self.beta2 = 0.01;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function r = getRange( self )
      r = [0,1];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function y = eval( self, alpha )
      t1 = sqrt(1+self.beta1)-self.beta1;
      t2 = sqrt(1+self.beta2)-self.beta2;
      y  = t1*sqrt( (1-alpha).^2 + self.beta2^2 ) + ...
           t2*sqrt( alpha.^2 + self.beta1^2 );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_D( self, alpha )
      t1 = sqrt(1+self.beta1)-self.beta1;
      t2 = sqrt(1+self.beta2)-self.beta2;
      Dy = t1*(alpha-1)./sqrt( (1-alpha).^2 + self.beta2^2 ) + ...
           t2*alpha./sqrt( alpha.^2 + self.beta1^2 );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = eval_DD( self, alpha )
      t1  = sqrt(1+self.beta1)-self.beta1;
      t2  = sqrt(1+self.beta2)-self.beta2;
      DDy = t1*self.beta2^2./( (1-alpha).^2 + self.beta2^2 )^(3/2) + ...
            t2*self.beta1^2./( alpha.^2 + self.beta1^2 )^(3/2);
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
