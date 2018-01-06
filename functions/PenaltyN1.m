classdef PenaltyN1 < FunctionND
  %
  % The Penalty Function #1
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

    function self = PenaltyN1( varargin )
      if nargin == 0
        n = int32(2) ;
      elseif nargin == 1
        n = varargin{1} ;          
      else
        error('PenaltyN1: too much argument in constructor') ;
      end
      if ~isinteger(n)
        error('PenaltyN1: argument must be an integer, found %s',class(n));
      end
      if n <= 1
        error('PenaltyN1: argument must be an integer > 1, found %d',n);
      end
      self@FunctionND(int32(n)) ;
      self.exact_solutions = zeros(n,0) ; % unknown solution 
      self.guesses         = (1:n).';
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      ap = 0.00001;
      t1 = dot( x, x ) - 0.25 ;
      t2 = sum ( ( x - 1.0 ).^2 );
      f  = ap * t2 + t1^2;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      ap = 0.00001;
      t1 = dot( x, x ) - 0.25 ;
      g  = (2.0 * ap) * ( x' - 1.0 ) + (4.0 * t1) * x' ;
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      ap = 0.00001;

      t1 = - 0.25 + dot( x, x );

      d1 = 2.0 * ap;
      th = 4.0 * t1;

      for i = 1 : self.N
        h(i,i) = d1 + th + 8.0 * x(i)^2;
        for j = 1 : i - 1
          h(i,j) = 8.0 * x(i) * x(j);
        end
      end

      for i = 1 : self.N
        for j = i + 1 : self.N
          h(i,j) = h(j,i);
        end
      end
    end
  end
end
