classdef Leon < FunctionND
  %
  % The Leon cubic valley function.
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

  methods

    function self = Leon()
      self@FunctionND(int32(2)) ;
      self.exact_solutions = [ 1.0; 1.0  ];     % one known solution 
      self.guesses         = [ -1.2; -1.0 ];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f1 = x(2) - x(1) * x(1) * x(1);
      f2 = 1.0 - x(1);
      f  = 100.0 * f1 * f1 + f2 * f2;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g    = zeros ( 1, 2 );
      g(1) = - 600.0 * ( x(2) - x(1)^3 ) * x(1) * x(1) - 2.0 * ( 1.0 - x(1) );
      g(2) = 200.0 * ( x(2) - x(1)^3 );
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( 2, 2 );
      h(1,1) = - 1200.0 * x(1) * x(2) + 3000.0 * x(1)^4 + 2.0;
      h(1,2) = - 600.0 * x(1) * x(1);
      h(2,1) = - 600.0 * x(1) * x(1);
      h(2,2) = 200.0;
    end
  end
end
