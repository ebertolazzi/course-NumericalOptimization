classdef Powell3D < FunctionND
  %
  % The Bohachevsky Function #1. (function for problem 38)
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
  %   Richard Brent,
  %   Algorithms for Minimization with Derivatives,
  %   Dover, 2002,
  %   ISBN: 0-486-41998-3,
  %   LC: QA402.5.B74.
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

    function self = Powell3D()
      self@FunctionND(int32(3)) ;
      self.exact_solutions = ones(3,1) ; % one known solution 
      self.guesses         = [ 0.0; 1.0; 2.0 ];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      if ( x(2) == 0.0 )
        term = 0.0;
      else
        arg  = ( x(1) + 2.0 * x(2) + x(3) ) / x(2);
        term = exp ( - arg^2 );
      end

      f = 3.0 ...
        - 1.0 / ( 1.0 + ( x(1) - x(2) )^2 ) ...
        - sin ( 0.5 * pi * x(2) * x(3) ) ...
        - term;
    end

    function g = grad( self, X )
      % use analitic gradient
      self.check_x(X);
      g = zeros ( 1, 3 );
      x = X(1) ;
      y = X(2) ;
      z = X(3) ;

      t1   = (x - y);
      t2   = (t1 ^ 2);
      t4   = (1 + t2) ^ 2;
      t5   = 1 / t4;
      t11  = cos(pi * y * z / 0.2e1);
      t12  = t11 * pi;
      g(1) = 2 * t5 * t1;
      g(2) = (-(2 * t5 * t1) - t12 * z / 0.2e1);
      g(3) = -(t12 * y / 0.2e1);

      if ( y ~= 0.0 )
        t1 = 0.1e1 / y;
        t5 = exp((-x - 0.20e1 * y - z) * t1);
        t6 = t5 * t1;
        t9 = y ^ 2;
        g(1) = g(1) - t6 ;
        g(2) = g(2) + (x + z) * t5 / t9;
        g(3) = g(3)-t6;
      end
    end

    function h = hessian( self, X )
      % use analitic hessian
      self.check_x(X);
      h = zeros ( 3, 3 );
      x = X(1) ;
      y = X(2) ;
      z = X(3) ;

      t2 = ((x - y) ^ 2);
      t3 = (x ^ 2);
      t4 = (x * y);
      t6 = (y ^ 2);
      t7 = 1 + t3 - 2 * t4 + t6;
      t8 = t7 ^ 2;
      t10 = 1 / t8 / t7;
      t12 = 8 * t2 * t10;
      t14 = (1 + t2) ^ 2;
      t16 = 2 / t14;
      t22 = (6 * t3 - 12 * t4 + 6 * t6 - 2) * t10;
      t25 = pi * y * z / 0.2e1;
      t26 = sin(t25);
      t27 = pi ^ 2;
      t28 = t26 * t27;
      t29 = z ^ 2;
      t36 = cos(t25);
      t40 = pi * (t26 * pi * y * z - 0.2e1 * t36) / 0.4e1;
      h(1,1) = -t12 + t16;
      h(1,2) = t22;
      h(1,3) = 0;
      h(2,1) = t22;
      h(2,2) = (-t12 + t16 + t28 * t29 / 0.4e1);
      h(2,3) = t40;
      h(3,1) = 0;
      h(3,2) = t40;
      h(3,3) = (t28 * t6 / 0.4e1);

      if ( y ~= 0.0 )
        t1 = y ^ 2;
        t7 = exp((-x - 0.20e1 * y - z) / y);
        t8 = 0.1e1 / t1 * t7;
        t13 = t7 * (-y + x + z) / t1 / y;
        t18 = t1 ^ 2;
        h(1,1) = t8;
        h(1,2) = -t13;
        h(1,3) = t8;
        h(2,1) = -t13;
        h(2,2) = (x + z) * t7 * (-0.2e1 * y + x + z) / t18;
        h(2,3) = -t13;
        h(3,1) = t8;
        h(3,2) = -t13;
        h(3,3) = t8;
      end
    end
  end
end
