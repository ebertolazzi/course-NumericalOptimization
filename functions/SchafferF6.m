classdef SchafferF6 < FunctionND

  %
  % The Schaffer Function F6.
  %
  % Reference:
  % @book{brent1973algorithms,
  %   title={Algorithms for Minimization Without Derivatives},
  %   author={Brent, Richard P.},
  %   isbn={9780486419985},
  %   lccn={01047459},
  %   series={Dover Books on Mathematics},
  %   year={1973},
  %   publisher={Dover Publications}
  % }
  %
  % Author:
  %   Enrico Bertolazzi
  %   Dipartimento di Ingegneria Indutriale
  %   Universita` degli Studi di Trento
  %   email: enrico.bertolazzi@unitn.it
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %

  methods (Hidden = true)

    function res = Sinc( ~, x )
      % function $ \sin(x)/x $ with taylor expansion 
      if abs(x) < 0.002
        res = 1+x*(1/6-x*x/20);
      else
        res = sin(x)/x;
      end
    end
  
    function res = Cosc( ~, x )
      % function $ 1-\cos(x)/x $ with taylor expansion 
      x2 = x*x ;
      if abs(x) < 0.02
        res = (1.0/2.0+(-.01/24.0+(1.0/720.0)*x2)*x2)*x;
      else
        res = (1-cos(x))/x;
      end
    end

    function res = dSinc( self, x )
      % function $ (1-\mathrm{sinc}(x))/x^2 $ with taylor expansion 
      if abs(x) < 0.018
        x2  = x*x ;
        res = (1-(x2/20)*(1-x2/42))/6 ;
      else
        res = (1-self.Sinc(x))/x^2 ;
      end
    end

    function res = dCosc( self, x )
      % function $ (x/2-\mathrm{cosc}(x))/x^3 $ with taylor expansion 
      if abs(x) < 0.026
        res = (1-(x2/30)*(1+x2/56))/24 ;
      else
        res = (1/2-self.Cosc(x)/x)/x^2 ;
      end
    end

    function varargout = a2( ~, r2 )
      a = 1/(1 + 0.001*r2);
      if nargout > 0
        varargout{1} = a*a ;
      end
      if nargout > 1
        varargout{2} = -0.002*varargout{1}*a ;
      end
      if nargout > 2
        varargout{3} = -0.003*varargout{2}*a ;
      end
    end

    function varargout = b( self, r2 )
      r = sqrt(r2);
      if nargout > 0
        varargout{1} = sin(r)^2-0.5 ;
      end
      if nargout > 1
        varargout{2} = self.Sinc(r)*cos(r) ;
      end
      if nargout > 2
        DS = self.dSinc(r) ;
        DC = self.dCosc(r) ;
        t1 = (2*DS-3)/4 ;
        t2 = (4*DC+6*DS+1)/8 ;
        t3 = ( DS*(DC-DS) - DC ) / 2 ;
        t4 = (DC^2)/2;
        varargout{3} = t1+(t2+(t3+t4*r2)*r2)*r2;
      end
    end

  end

  methods

    function self = SchafferF6()
      self@FunctionND(int32(2)) ;
      self.exact_solutions = [ 0 ; 0 ];       % one known solution 
      self.guesses         = [ -5.0 ; +10.0 ] ; % one guess
    end

    function f = eval(self,x)
      % evaluate Rosenbrock (2D) function.
      self.check_x(x);
      X  = x(1) ;
      Y  = x(2) ;
      r2 = X^2+Y^2;
      a2 = self.a2(r2);
      b  = self.b(r2);
      f  = 0.5 + a2*b;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      X         = x(1) ;
      Y         = x(2) ;
      r2        = X^2+Y^2;
      [a2,a2_1] = self.a2(r2);
      [b,b_1]   = self.b(r2);
      d1        = 2*(a2*b_1+a2_1*b) ;
      g         = d1 * [X,Y] ;
    end

    function h = hessian( self, x )
      % use analitic hessian
      h = zeros( 2, 2 ) ;
      self.check_x(x);
      X              = x(1) ;
      Y              = x(2) ;
      X2             = X^2;
      Y2             = Y^2;
      r2             = X2+Y2;
      [a2,a2_1,a2_2] = self.a2(r2);
      [b,b_1,b_2]    = self.b(r2);
      d1             = 2*(a2*b_1+a2_1*b) ;
      d2             = 4*(a2*b_2+2*a2_1*b_1+a2_2*b) ;

      h(1,1) = (d2*X2+d1)*r2;
      h(2,2) = (d2*Y2+d1)*r2;
      h(2,1) = d2*r2*X*Y;
      h(1,2) = h(2,1);
    end
  end
end
