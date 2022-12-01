classdef Hilbert < FunctionND
  %
  % The Hilbert Matrix Function F = x'Ax
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
    function self = Hilbert( varargin )
      if nargin == 0
        n = int32(2);
      elseif nargin == 1
        n = varargin{1};
      else
        error('Hilbert: too much argument in constructor');
      end
      if ~isinteger(n)
        error('Hilbert: argument must be an integer, found %s',class(n));
      end
      if n <= 0
        error('Hilbert: argument must be an integer > 0, found %d',n);
      end
      self@FunctionND(int32(n));
      self.exact_solutions = zeros( n, 1 );
      self.guesses = ones( n, 1 );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = 0.0;
      for i = 1 : self.N
        for j = 1 : self.N
          f = f + x(i) * x(j) / double( i + j - 1 );
        end
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, self.N );
      for i = 1 : self.N
        g(i) = 0.0;
        for j = 1 : self.N
          g(i) = g(i) + 2.0 * x(j) / double( i + j - 1 );
        end
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros( self.N, self.N );
      for i = 1 : self.N
        for j = 1 : self.N
          h(i,j) = 2.0 / double( i + j - 1 );
        end
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
