classdef GoldenSearch < LinesearchArmijo
  properties (SetAccess = private, Hidden = true)
    tau
    tol      % tolerance for |b-a|
    max_iter % maximum number of admitted iterate for Golde search
  end
  methods
    function self = GoldenSearch()
      % PASSA al costruttore della super-classe
      self@LinesearchArmijo() ;
      self.tau      = (sqrt(5)-1)/2 ;
      self.tol      = 1e-3 ;
      self.max_iter = 10 ;
    end

    function [a,b] = minimize( self, a_in, b_in )
      % controllo che b > al
      a      = a_in ;
      b      = b_in ;
      dlen   = self.tau*(b-a) ;
      lambda = b-dlen ;
      mu     = a+dlen ;
      fl     = self.fun1D.eval(lambda);
      fm     = self.fun1D.eval(mu);
      tolerance = self.tol*(b-a); % tol viene dalla classe Minimization1D
      while (b-a) > tolerance
        if fl > fm
          % seleziono intervallo [lambda,b]
          a      = lambda ;
          lambda = mu;      fl = fm ;
          dlen   = self.tau*(b-a) ;
          mu     = a+dlen ;
          fm     = self.fun1D.eval(mu);
        else
          % seleziono intervallo [a,mu]
          b      = mu;
          mu     = lambda; fm = fl ;
          dlen   = self.tau*(b-a) ;
          lambda = b-dlen ;
          fl     = self.fun1D.eval(lambda);
        end
      end
    end
  end
end
