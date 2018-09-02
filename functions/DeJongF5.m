classdef DeJongF5 < FunctionND
  %
  % The De Jong Function F5
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

    function self = DeJongF5()
      self@FunctionND(int32(2));
      self.exact_solutions = [ -32.0; -32.0 ];  % one known solution
      self.guesses         = [ -32.01; -32.02 ];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      jroot = 5;
      k = 500.0;

      fi = k;

      for j = 1 : jroot * jroot

        j1 = mod ( j - 1, jroot ) + 1;
        a1 = -32 + j1 * 16;

        j2 = ( j - 1 ) / jroot;
        a2 = -32 + j2 * 16;

        fj = j + ( x(1) - a1 )^6 + ( x(2) - a2 )^6;

        fi = fi + 1.0 / fj;

      end

      f = 1.0 / fi;
    end


    function g = grad( self, x )
      g = self.FD_grad( x );
    end

    function h = hessian( self, x )
      g = self.FD_hessianj( x );
    end

  end
end
