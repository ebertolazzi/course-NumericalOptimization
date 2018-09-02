classdef BraninRCOS < FunctionND
  %
  % The Branin RCOS Function
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
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

    function self = BraninRCOS()
      self@FunctionND(int32(2));
      self.exact_solutions = [ -pi, 12.275;
                                pi,  2.275;
                                9.42478, 2.475 ].'; % 3 known solutions
      self.guesses = [ -1.0; 1.0 ];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      a = 1.0;
      d = 6.0;
      e = 10.0;

      b = 5.1 / ( 4.0 * pi^2 );
      c = 5.0 / pi;
      ff = 1.0 / ( 8.0 * pi );

      f = a * ( x(2) - b * x(1)^2 + c * x(1) - d )^2 + e * ( 1.0 - ff ) * cos ( x(1) ) + e;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 2 );

      a = 1.0;
      d = 6.0;
      e = 10.0;

      b = 5.1 / ( 4.0 * pi^2 );;
      c = 5.0 / pi;
      ff = 1.0 / ( 8.0 * pi );

      tmp = 2.0 * a * ( x(2) - b * x(1)^2 + c * x(1) - d );

      g(1) = tmp * ( c - 2 * b * x(1) ) - e * ( 1.0 - ff ) * sin ( x(1) );
      g(2) = tmp;
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h  = zeros ( 2, 2 );
      a = 1.0;
      d = 6.0;
      e = 10.0;

      b  = 5.1 / ( 4.0 * pi^2 );
      c  = 5.0 / pi;
      ff = 1.0 / ( 8.0 * pi );

      tmp = c - 2.0 * b * x(1);

      h(1,1) = 2.0 * a * tmp^2 - 4.0 * a * b * ( x(2) - b * x(1)^2 + c * x(1) - d ) - e * ( 1.0 - ff ) * cos ( x(1) );
      h(1,2) = 2.0 * a * tmp;
      h(2,1) = h(1,2);
      h(2,2) = 2.0 * a;
    end
  end
end
