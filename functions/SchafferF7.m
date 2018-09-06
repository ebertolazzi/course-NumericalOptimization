classdef SchafferF7 < FunctionND

  %
  % The Schaffer Function F7.
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

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = SchafferF7()
      self@FunctionND(int32(2));
      self.exact_solutions = [ 0; 0 ];       % one known solution
      self.guesses         = [ -5.0; +10.0 ]; % one guess
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      % evaluate Rosenbrock (2D) function.
      self.check_x(x);
      X = x(1);
      Y = x(2);
      r = hypot( X, Y );
      f = sqrt ( r ) * ( 1.0 + ( sin ( 50.0 * r^0.2 ) )^2 );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      X = x(1);
      Y = x(2);
      r = hypot( X, Y );
      g = zeros(1,2);

      if r > 0
        a  = sqrt ( r );
        ar = 0.5 / sqrt ( r );

        b  = 1.0 + ( sin ( 50.0 * r^0.2 ) )^2;
        br = 10.0 * sin ( 100.0 * r^0.2 ) * r^(-0.8);

        rx1 = X / r;
        rx2 = Y / r;

        g(1) = ( ar * b + a * br ) * rx1;
        g(2) = ( ar * b + a * br ) * rx2;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, x )
      % use analitic hessian
      % use analitic gradient
      self.check_x(x);
      X = x(1);
      Y = x(2);

      h = zeros ( 2, 2 );
      r = hypot( X, Y );

      rx1 = X / r;
      rx2 = Y / r;

      r3    = r^3;
      rx1x1 = Y^2 / r3;
      rx1x2 = - X*Y / r3;
      rx2x2 = X^2 / r3;
      %
      %  F = A * B
      %  dFdX1 = ( Ar * B + A * Br ) * Rx1
      %  d2FdX1dX1 = ( Arr * B + Ar * Br ) * Rx1^2 + ( Ar * B + A * Br ) * Rx1x1
      %  etc
      %
      a   = sqrt ( r );
      ar  = 0.5 / a;
      arr = - 0.25 / (r*a);

      S50  = sin( 50.0 * r^0.2 );
      S100 = sin( 100.0 * r^0.2 );

      b   = 1 + S50^2;
      br  = 10 * S100 * r^(-0.8);
      brr = 200 * cos ( 100 * r^0.2 ) * r^(-1.6) - 8.0 * S100 * r^(-1.8);

      t1 = arr * b + 2 * ar * br + a * brr;
      t2 = ar * b + a * br;

      h(1,1) = t1 * rx1 * rx1 + t2 * rx1x1;
      h(1,2) = t1 * rx1 * rx2 + t2 * rx1x2;
      h(2,1) = h(1,2);
      h(2,2) = t1 * rx2 * rx2 + t2 * rx2x2;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = eval_FG( self, x )
      f = self.eval(x);
      g = self.grad(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g,H] = eval_FGH( self, x )
      f = self.eval(x);
      g = self.grad(x);
      H = self.hessian(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
