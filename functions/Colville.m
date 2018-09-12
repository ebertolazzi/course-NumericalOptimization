classdef Colville < FunctionND
  %
  % The Colville Polynomial
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
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Colville()
      self@FunctionND(int32(4));
      self.exact_solutions = [ 1; 1; 1; 1 ];     % one known solution
      self.guesses         = [ 0.5; 1.0; -0.5; -1.0 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = 100.0 * ( x(2) - x(1)^2 )^2 ...
        + ( 1.0 - x(1) )^2 ...
        + 90.0 * ( x(4) - x(3)^2 )^2 ...
        + ( 1.0 - x(3) )^2 ...
        + 10.1 * ( ( x(2) - 1.0 )^2 + ( x(4) - 1.0 )^2 ) ...
        + 19.8 * ( x(2) - 1.0 ) * ( x(4) - 1.0 );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 4 );

      g(1) = 400.0 * x(1)^3 - 400.0 * x(2) * x(1) + 2.0 * x(1) - 2.0;
      g(2) = -200.0 * x(1)^2 + 220.2 * x(2) + 19.8 * x(4) - 40.0;
      g(3) = -360.0 * x(3) * x(4) + 360.0 * x(3)^3 + 2.0 * x(3) - 2.0;
      g(4) = + 180.0 * x(4) - 180.0 * x(3)^2 + 20.2 * x(4) + 19.8 * x(2) - 40.0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( 4, 4 );
      h(1,1) = 1200.0 * x(1)^2 - 400.0 * x(2) + 2.0;
      h(1,2) = - 400.0 * x(1);

      h(2,1) = -400.0 * x(1);
      h(2,2) = 220.2;
      h(2,4) = 19.8;

      h(3,3) = -360.0 * x(4) + 1080.0 * x(3)^2 + 2.0;
      h(3,4) = - 360.0 * x(3);

      h(4,2) = 19.8;
      h(4,3) = - 360.0 * x(3);
      h(4,4) = 200.2;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = eval_FG( self, x )
      f = self.eval(x);
      g = self.grad(x);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g,H] = eval_FGH( self, x )
      f = self.eval(x);
      g = self.grad(x);
      H = self.hessian(x);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
