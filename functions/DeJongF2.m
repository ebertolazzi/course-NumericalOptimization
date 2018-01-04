classdef DeJongF2 < FunctionND
  %
  % The De Jong Function F2
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

    function self = DeJongF2()
      self@FunctionND(int32(2)) ;
      self.exact_solutions = [ 1 ; 1 ];     % one known solution 
      self.guesses         = [ -2.048 ; 2.048 ];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = 100.0 * ( x(1) * x(1) - x(2) )^2 + ( 1.0 - x(1) )^2;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g    = zeros ( 1, 2 );
      g(1) = 400.0 * ( x(1) * x(1) - x(2) ) * x(1) - 2.0 + 2.0 * x(1);
      g(2) = - 200.0 * ( x(1) * x(1) - x(2) );
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( 2, 2 );
      h(1,1) = 1200.0 * x(1) * x(1) - 400.0 * x(2) + 2.0;
      h(1,2) = - 400.0 * x(1);
      h(2,1) = -400.0 * x(1);
      h(2,2) = 200.0;
    end
  end
end
