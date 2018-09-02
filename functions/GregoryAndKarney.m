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

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = x(1) * x(1) + 2.0 * dot( x, x );
      for i = 1 : self.N - 1
        f = f - 2.0 * x(i) * x(i+1);
      end
      f = f - 2.0 * x(1);
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, self.N );

      for i = 1 : self.N

        if ( i == 1 )
          g(i) = x(i);
        else
          g(i) = 2.0 * x(i);
        end

        if ( 1 < i )
          g(i) = g(i) - x(i-1);
        end

        if ( i < self.N )
          g(i) = g(i) - x(i+1);
        end

      end

      g(1) = g(1) - 2.0;
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( self.N, self.N );

      for i = 1 : self.N
        for j = 1 : self.N
          if i == j
            if i == 1
              h(i,j) = 2.0;
            else
              h(i,j) = 4.0;
            end
          elseif i == j + 1 || i == j - 1
            h(i,j) = - 2.0;
          else
            h(i,j) = 0.0;
          end
        end
      end
    end
  end
end
