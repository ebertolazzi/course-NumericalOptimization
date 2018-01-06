classdef Box3 < FunctionND
  %
  % The Box 3-dimensional function.
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

    function self = Box3()
      self@FunctionND(int32(3)) ;
      self.exact_solutions = [ 1.0; 10.0; 1.0];     % one known solution 
      self.guesses         = [ 0.0; 10.0; 5.0 ];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      f = 0.0;
      for i = 1 : 10
        c = - i / 10.0;
        fi = exp ( c * x(1) ) - exp ( c * x(2) ) - x(3) * ( exp ( c ) - exp ( 10.0 * c ) );
        f = f + fi * fi;
      end
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 3 );

      for i = 1 : 10
        c = - i / 10.0;
        fi = exp ( c * x(1) ) - exp ( c * x(2) ) - x(3) * ( exp ( c ) - exp ( 10.0 * c ) );

        dfidx1 =   c * exp ( c * x(1) );
        dfidx2 = - c * exp ( c * x(2) );
        dfidx3 = - ( exp ( c ) - exp ( 10.0 * c ) );

        g(1) = g(1) + 2.0 * fi * dfidx1;
        g(2) = g(2) + 2.0 * fi * dfidx2;
        g(3) = g(3) + 2.0 * fi * dfidx3;

      end
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( 3, 3 );

      for i = 1 : 10

        c = - i / 10.0;

        fi = exp ( c * x(1) ) - exp ( c * x(2) ) - x(3) * ( exp ( c ) - exp ( 10.0 * c ) );;

        dfidx1   =   c     * exp ( c * x(1) );
        d2fidx11 =   c * c * exp ( c * x(1) );
        dfidx2   = - c     * exp ( c * x(2) );
        d2fidx22 = - c * c * exp ( c * x(2) );
        dfidx3   = - ( exp ( c ) - exp ( 10.0 * c ) );

        h(1,1) = h(1,1) + 2.0 * dfidx1 * dfidx1 + 2.0 * fi * d2fidx11;
        h(1,2) = h(1,2) + 2.0 * dfidx1 * dfidx2;
        h(1,3) = h(1,3) + 2.0 * dfidx1 * dfidx3;

        h(2,1) = h(2,1) + 2.0 * dfidx2 * dfidx1;
        h(2,2) = h(2,2) + 2.0 * dfidx2 * dfidx2 + 2.0 * fi * d2fidx22;
        h(2,3) = h(2,3) + 2.0 * dfidx2 * dfidx3;

        h(3,1) = h(3,1) + 2.0 * dfidx3 * dfidx1;
        h(3,2) = h(3,2) + 2.0 * dfidx3 * dfidx2;
        h(3,3) = h(3,3) + 2.0 * dfidx3 * dfidx3;

      end
    end
  end
end
