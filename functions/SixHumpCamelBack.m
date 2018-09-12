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
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = SixHumpCamelBack()
      self@FunctionND(int32(2));
      self.exact_solutions = [ -0.0898,  0.7126;
                                0.0898, -0.7126 ].';
      self.guesses = [ -1.5; 0.5 ]; % one guess
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval( self, xx )
      self.check_x( xx );
      x  = xx(1);
      y  = xx(2);
      x2 = x*x;
      y2 = y*y;
      f  = ( 4 + x2*(x2/3-2.1) ) * x2 + x * y + 4*( y2 - 1 ) * y2;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, xx )
      % use analitic gradient
      self.check_x( xx );
      x  = xx(1);
      y  = xx(2);
      x2 = x*x;
      g  = zeros( 1, 2 );
      g(1) = 2*x*(4+x2*(x2-4.2)) + y;
      g(2) = 2*y*(8*y*y-4) + x;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, xx )
      % use analitic hessian
      self.check_x( xx );
      x  = xx(1);
      y  = xx(2);
      x2 = x*x;
      y2 = y*y;

      h = zeros( 2, 2 );
      h(1,1) = (10*x2-25.2)*x2+8;
      h(1,2) = 1;
      h(2,1) = 1;
      h(2,2) = 48*y2-8;
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
