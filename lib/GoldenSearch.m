classdef GoldenSearch < LinesearchForwardBackward
  
    properties (SetAccess = private, Hidden = true)
    tau
    tol      % tolerance for |b-a|
    max_iter % maximum number of admitted iterate for Golde search
  end
  
  methods
  
    function self = GoldenSearch()
      % PASSA al costruttore della super-classe
      self@LinesearchForwardBackward() ;
      self.tau      = (sqrt(5)-1)/2 ;
      self.tol      = 1e-3 ;
      self.max_iter = 10 ;
    end

    function setMaxIteration( self, max_iter )
      if length(max_iter) > 1 || ~isinteger(max_iter)
        error('GoldenSearch, expected a scalar  integer\n') ;
      end
      if max_iter < 0 || max_iter > 10000
        error('GoldenSearch, bad number of iterator %d\n',max_iter) ;
      end
      self.max_iter = max_iter ;
    end

    function setTolerance( self, tol )
      if tol <= 0
        error('GoldenSearch, bad tolerance %g\n',tol) ;
      end
      self.tol = tol ;
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
      for k=1:self.max_iter
        if (b-a) < tolerance ; break ; end
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

    function [alpha,ok] = search( self, alpha_guess )
      [~,alpha1,ierr] = self.ForwardBackward( alpha_guess );
      ok = ierr >= 0 ;
      if ok
        [a,b] = self.minimize( 0, alpha1 ) ;
        alpha = (a+b)/2 ;
      else
        alpha = alpha_guess/100 ;
      end    
    end
    
  end
end
