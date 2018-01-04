classdef Gauss < FunctionND
  %
  % The Gaussian function.
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
  %   Zbigniew Michalewicz,
  %   Genetic Algorithms + Data Structures = Evolution Programs,
  %   Third Edition,
  %   Springer Verlag, 1996,
  %   ISBN: 3-540-60676-9,
  %   LC: QA76.618.M53.
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

  properties (SetAccess = protected, Hidden = true)
    y
  end

  methods

    function self = Gauss()
      self@FunctionND(int32(3)) ;
      self.exact_solutions = [ 0 ; 0 ; 0 ];     % one known solution 
      self.guesses         = [ 0.4; 1.0; 0.0 ];

      self.y = [ 0.0009, 0.0044, 0.0175, 0.0540, ...
                 0.1295, 0.2420, 0.3521, 0.3989, ...
                 0.3521, 0.2420, 0.1295, 0.0540, ...
                 0.0175, 0.0044, 0.0009 ] ;
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = 0.0;
      for i = 1 : 15
        t = x(1) * exp ( - 0.5 * x(2) * ( 3.5 - 0.5 * ( i - 1 ) - x(3) )^2 ) - self.y(i);
        f = f + t * t;
      end
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 3 );
      for i = 1 : 15
        d1  = 0.5 * ( i - 1 );
        d2  = 3.5 - d1 - x(3);
        arg = - 0.5 * x(2) * d2 * d2;
        t   = x(1) * exp ( arg ) - self.y(i);

        g(1) = g(1) + 2.0 * exp ( arg ) * t;
        g(2) = g(2) - x(1) * exp ( arg ) * t * d2 * d2;
        g(3) = g(3) + 2.0 * x(1) * x(2) * exp ( arg ) * t * d2;
      end
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);

      h = zeros ( 3, 3 );

      for i = 1 : 15

        d1 = 0.5 * ( i - 1 );
        d2 = 3.5 - d1 - x(3);
        arg = 0.5 * x(2) * d2 * d2;
        r = exp ( - arg );
        t = x(1) * r - self.y(i);
        t1 = 2.0 * x(1) * r - self.y(i);

        h(1,1) = h(1,1) + r * r;
        h(2,2) = h(2,2) + r * t1 * d2^4;
        h(3,3) = h(3,3) + r * ( x(2) * t1 * d2 * d2 - t );
        h(2,1) = h(2,1) - r * t1 * d2 * d2;
        h(3,1) = h(3,1) + d2 * r * t1;
        h(3,2) = h(3,2) + d2 * r * ( t - arg * t1 );

      end

      h(1,1) = 2.0 * h(1,1);
      h(2,2) = 0.5 * x(1) * h(2,2);
      h(3,3) = 2.0 * x(1) * x(2) * h(3,3);
      h(3,1) = 2.0 * x(2) * h(3,1);
      h(3,2) = 2.0 * x(1) * h(3,2);

      h(1,2) = h(2,1);
      h(1,3) = h(3,1);
      h(2,3) = h(3,2);
    end
  end
end
