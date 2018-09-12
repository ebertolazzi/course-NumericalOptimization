classdef GregoryAndKarney < FunctionND
  %
  % The Gregory and Karney Tridiagonal Matrix Function.
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
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
    function self = GregoryAndKarney( varargin )
      if nargin == 0
        n = int32(2);
      elseif nargin == 1
        n = varargin{1};
      else
        error('GregoryAndKarney: too much argument in constructor');
      end
      if ~isinteger(n)
        error('GregoryAndKarney: argument must be an integer, found %s',class(n));
      end
      if n <= 0
        error('GregoryAndKarney: argument must be an integer > 0, found %d',n);
      end
      self@FunctionND(int32(n));

      self.exact_solutions = zeros( n, 1 );
      for i = 1:n
        self.exact_solutions(i) = n + 1 - i;
      end
      self.guesses = zeros( n, 1 );

    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval( self, x )
      % evaluate function
      self.check_x(x);
      f = x(1) * x(1) + 2 * dot( x, x );
      for i = 1 : self.N - 1
        f = f - 2.0 * x(i) * x(i+1);
      end
      f = f - 2.0 * x(1);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros( 1, self.N );
      g(1) = 6*x(1);
      for i=2:self.N
        g(i) = 4*x(i);
      end
      for i = 1 : self.N - 1
        g(i)   = g(i)   - 2*x(i+1);
        g(i+1) = g(i+1) - 2*x(i);
      end
      g(1) = g(1) - 2;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros( self.N, self.N );
      h(1,1) = 6;
      for i = 1:self.N
        h(i,i) = 4;
      end
      for i = 1:self.N-1
        h(i,i+1) = h(i,i+1)-2;
        h(i+1,i) = h(i+1,i)-2;
      end
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
