classdef Himmelblau < FunctionND
  %
  % The Himmelblau function
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
  %   David Himmelblau,
  %   Applied Nonlinear Programming,
  %   McGraw Hill, 1972,
  %   ISBN13: 978-0070289215,
  %   LC: T57.8.H55.
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

    function self = Himmelblau()
      self@FunctionND(int32(2)) ;
      self.exact_solutions = [ 3.0; 2.0];     % one known solution 
      self.guesses         = [ -1.3; 2.7];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = ( x(1)^2 + x(2) - 11.0 )^2 + ( x(1) + x(2)^2 - 7.0 )^2;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 2 );
      g(1) = 4.0 * ( x(1)^2 + x(2) - 11.0 ) * x(1) + 2.0 * ( x(1) + x(2)^2 - 7.0 );
      g(2) = 2.0 * ( x(1)^2 + x(2) - 11.0 ) + 4.0 * ( x(1) + x(2)^2 - 7.0 ) * x(2);
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( 2, 2 ) ;
      h(1,1) = 8.0 * x(1) * x(1) + 4.0 * ( x(1) * x(1) + x(2) - 11.0 ) + 2.0;
      h(1,2) = 4.0 * x(1) + 4.0 * x(2) ;
      h(2,1) = 4.0 * x(1) + 4.0 * x(2) ;
      h(2,2) = 2.0 + 8.0 * x(2) * x(2) + 4.0 * ( x(1) + x(2) * x(2) - 7.0 );
    end
  end
end
