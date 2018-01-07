classdef SixHumpCamelBack < FunctionND

  %
  % The Schaffer Function F6.
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

    function self = SixHumpCamelBack()
      self@FunctionND(int32(2)) ;
      self.exact_solutions = [ -0.0898,  0.7126 ;
                                0.0898, -0.7126 ].'; 
      self.guesses = [ -1.5; 0.5 ] ; % one guess
    end

    function f = eval(self,x)
      % evaluate Rosenbrock (2D) function.
      self.check_x(x);
      f = ( 4.0 - 2.1 * x(1)^2 + x(1)^4 / 3.0 ) * x(1)^2 + x(1) * x(2) + 4.0 * ( x(2)^2 - 1.0 ) * x(2)^2;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 2 );
      g(1) = 8.0 * x(1) - 8.4 * x(1)^3 + 2.0 * x(1)^5 + x(2);
      g(2) = x(1) - 8.0 * x(2) + 16.0 * x(2)^3;
    end

    function h = hessian( self, x )
      % use analitic hessian
      h = zeros ( 2, 2 );
      h(1,1) = 8.0 - 25.2 * x(1)^2 + 10.0 * x(1)^4;
      h(1,2) = 1.0;
      h(2,1) = 1.0;
      h(2,2) = -8.0 + 48.0 * x(2)^2;
    end
  end
end
