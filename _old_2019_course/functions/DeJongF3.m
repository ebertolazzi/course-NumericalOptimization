classdef DeJongF3 < FunctionND
  %
  % The De Jong Function F3
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
    function self = DeJongF3()
      self@FunctionND(int32(5.12));
      self.exact_solutions = [ -5; -5; -5; -5; -5 ]; % one known solution
      self.guesses         = [ -5.12; -2.56; 0; 2.56; 5.12 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = sum( floor (x) );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, x )
      g = self.FD_grad( x );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function H = hessian( self, x )
      H = self.FD_hessianj( x );
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
