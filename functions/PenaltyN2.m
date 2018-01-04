classdef PenaltyN2 < FunctionND
  %
  % The Penalty Function #2
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

    function self = PenaltyN2(n)
      if ~isinteger(n)
        error('Hilbert: argument must be an integer, found %s',class(n));
      end
      if n <= 1
        error('Hilbert: argument must be an integer > 1, found %d',n);
      end
      self@FunctionND(int32(n)) ;
      self.exact_solutions = zeros(n,0) ; % unknown solution 
      self.guesses         = 0.5 * ones ( n, 1 ) ;
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      ap = 0.00001;

      t1 = -1.0;
      t2 = 0.0;
      t3 = 0.0;
      d2 = 1.0;
      s2 = 0.0;

      for j = 1 : self.N
        t1 = t1 + ( self.N - j + 1 ) * x(j)^2;
        s1 = exp ( x(j) / 10.0 );
        if ( 1 < j )
          s3 = s1 + s2 - d2 * ( exp ( 0.1 ) + 1.0 );
          t2 = t2 + s3 * s3;
          t3 = t3 + ( s1 - 1.0 / exp ( 0.1 ) )^2;
        end
        s2 = s1;
        d2 = d2 * exp ( 0.1 );
      end

      f = ap * ( t2 + t3 ) + t1 * t1 + ( x(1) - 0.2 )^2;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, self.N );

      ap = 0.00001;

      t1 = -1.0;
      for j = 1 : self.N
        t1 = t1 + ( self.N - j + 1 ) * x(j)^2;
      end

      d2 = 1.0;
      th = 4.0 * t1;
      s2 = 0.0;
      for j = 1 : self.N
        g(j) = ( self.N - j + 1 ) * x(j) * th;
        s1 = exp ( x(j) / 10.0 );
        if ( 1 < j )
          s3 = s1 + s2 - d2 * ( exp ( 0.1 ) + 1.0 );
          g(j) = g(j) + ap * s1 * ( s3 + s1 - 1.0 / exp ( 0.1 ) ) / 5.0;
          g(j-1) = g(j-1) + ap * s2 * s3 / 5.0;
        end
        s2 = s1;
        d2 = d2 * exp ( 0.1 );
      end

      g(1) = g(1) + 2.0 * ( x(1) - 0.2 );
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      ap = 0.00001;

      t1 = - 1.0;
      for j = 1 : self.N
        t1 = t1 + ( self.N - j + 1 ) * x(j) * x(j);
      end

      d1 = exp ( 0.1 );
      d2 = 1.0;
      s2 = 0.0;
      th = 4.0 * t1;
      
      h = zeros(self.N,self.N) ;

      for j = 1 : self.N

        h(j,j) = 8.0 * ( ( self.N - j + 1 ) * x(j) )^2 + ( self.N - j + 1 ) * th;

        s1 = exp ( x(j) / 10.0 );

        if ( 1 < j )

          s3 = s1 + s2 - d2 * ( d1 + 1.0 );
          h(j,j) = h(j,j) + ap * s1 * ( s3 + s1 - 1.0 / d1 + 2.0 * s1 ) / 50.0;
          h(j-1,j-1) = h(j-1,j-1) + ap * s2 * ( s2 + s3 ) / 50.0;
          for k = 1 : j - 1
            t1 = exp ( k / 10.0 );
            h(j,k) = 8.0 * ( n - j + 1 ) *  ( n - k + 1 ) * x(j) * x(k);
          end

          h(j,j-1) = h(j,j-1) + ap * s1 * s2 / 50.0;

        end

        s2 = s1;
        d2 = d1 * d2;

      end

      h(1,1) = h(1,1) + 2.0;

      for i = 1 : self.N
        h(i,i+1:n) = h(i+1:n,i);
      end
    end
  end
end
