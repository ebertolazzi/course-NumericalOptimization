classdef GoldenSearch < LinesearchArmijo
  properties (SetAccess = private, Hidden = true)
    tau
  end
  methods
    function self = GoldenSearch( varargin )
      % PASSA al costruttore della super-classe
      self@Minimization1D( varargin{:} ) ;
      self.tau = (sqrt(5)-1)/2 ;
    end

    function [a,b] = minimize( self, a_in, b_in )
      % controllo che b > al
      a      = a_in ;
      b      = b_in ;
      dlen   = self.tau*(b-a) ;
      lambda = b-dlen ;
      mu     = a+dlen ;
      fl     = self.funObj.eval(lambda);
      fm     = self.funObj.eval(mu);
      tolerance = self.tol*(b-a); % tol viene dalla classe Minimization1D
      while (b-a) > tolerance
        if fl > fm
          % seleziono intervallo [lambda,b]
          a      = lambda ;
          lambda = mu;      fl = fm ;
          dlen   = self.tau*(b-a) ;
          mu     = a+dlen ;
          fm     = self.funObj.eval(mu);
        else
          % seleziono intervallo [a,mu]
          b      = mu;
          mu     = lambda; fm = fl ;
          dlen   = self.tau*(b-a) ;
          lambda = b-dlen ;
          fl     = self.funObj.eval(lambda);
        end
      end
    end
  end
end
