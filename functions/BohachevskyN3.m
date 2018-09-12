classdef BohachevskyN3 < FunctionND
  %
  % The Bohachevsky Function #2. (function for problem 38)
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

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = BohachevskyN3()
      self@FunctionND(int32(2));
      self.exact_solutions = [ 0; 0 ];     % one known solution
      self.guesses         = [ 0.5; 1.0];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      a1 = 3.0 * pi * x(1);
      a2 = 4.0 * pi * x(2);
      f  = x(1)^2 + 2.0 * x(2)^2 - 0.3 * cos ( a1 + a2 ) + 0.3;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g  = zeros ( 1, 2 );
      a1 = 3.0 * pi * x(1);
      a2 = 4.0 * pi * x(2);
      g(1) = 2.0 * x(1) + 0.9 * pi * sin ( a1 + a2 );
      g(2) = 4.0 * x(2) + 1.2 * pi * sin ( a1 + a2 );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h  = zeros ( 2, 2 );
      a1 = 3.0 * pi * x(1);
      a2 = 4.0 * pi * x(2);
      h(1,1) = 2.0 + 2.7 * pi^2 * cos ( a1 + a2 );
      h(1,2) = 3.6 * pi^2 * cos ( a1 + a2 );
      h(2,1) = h(1,2);
      h(2,2) = 4.0 + 16.0*0.3* pi^2 * cos ( a1 + a2 );
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
