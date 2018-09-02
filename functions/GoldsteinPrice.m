classdef GoldsteinPrice < FunctionND
  %
  % The Goldstein Price Polynomial
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
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %

  methods

    function self = GoldsteinPrice()
      self@FunctionND(int32(2));
      self.exact_solutions = [ 0; -1 ];     % one known solution
      self.guesses         = [ -0.5; +0.25];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      a = x(1) + x(2) + 1.0;
      b = 19.0 - 14.0 * x(1) + 3.0 * x(1) * x(1) - 14.0 * x(2) + 6.0 * x(1) * x(2) + 3.0 * x(2) * x(2);
      c = 2.0 * x(1) - 3.0 * x(2);
      d = 18.0 - 32.0 * x(1) + 12.0 * x(1) * x(1) + 48.0 * x(2) - 36.0 * x(1) * x(2) + 27.0 * x(2) * x(2);
      f = ( 1.0 + a * a * b ) * ( 30.0 + c * c * d );
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 2 );

      a = x(1) + x(2) + 1.0;

      b = 19.0 - 14.0 * x(1) + 3.0 * x(1)^2 - 14.0 * x(2) + 6.0 * x(1) * x(2) + 3.0 * x(2)^2;

      dbdx1 = - 14.0 + 6.0 * x(1) + 6.0 * x(2);
      dbdx2 = - 14.0 + 6.0 * x(1) + 6.0 * x(2);

      c = 2.0 * x(1) - 3.0 * x(2);

      d = 18.0 - 32.0 * x(1) + 12.0 * x(1)^2 + 48.0 * x(2) - 36.0 * x(1) * x(2) + 27.0 * x(2)^2;
      dddx1 = - 32.0 + 24.0 * x(1) - 36.0 * x(2);
      dddx2 = 48.0 - 36.0 * x(1) + 54.0 * x(2);

      g(1) = ( 1.0 + a^2 * b ) * ( 4.0 * c * d + c^2 * dddx1 ) + ( 2.0 * a * b + a^2 * dbdx1 ) * ( 30.0 + c^2 * d );
      g(2) = ( 1.0 + a^2 * b ) * ( -6.0 * c * d + c^2 * dddx2 ) + ( 2.0 * a * b + a^2 * dbdx2 ) * ( 30.0 + c^2 * d );
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( 2, 2 );

      a = x(1) + x(2) + 1.0;

      b = 19.0 - 14.0 * x(1) + 3.0 * x(1)^2 - 14.0 * x(2) + 6.0 * x(1) * x(2) + 3.0 * x(2)^2;
      e = - 14.0 + 6.0 * x(1) + 6.0 * x(2);
      c = 2.0 * x(1) - 3.0 * x(2);
      d = 18.0 - 32.0 * x(1) + 12.0 * x(1)^2 + 48.0 * x(2) - 36.0 * x(1) * x(2) + 27.0 * x(2)^2;
      r = - 32.0 + 24.0 * x(1) - 36.0 * x(2);
      s = 48.0 - 36.0 * x(1) + 54.0 * x(2);

      h(1,1) = 2.0 * ( 2.0 * a * b + a^2 * e ) ...
        * ( 4.0 * c * d + c^2 * r ) + ( 1.0 + a^2 * b ) ...
        * ( 8.0 * d + 4.0 * c * r + 4.0 * c * r + 24.0 * c^2 ) ...
        + ( 2.0 * b + 2.0 * a * e + 2.0 * a * e + 6.0 * a^2 ) ...
        * ( 30.0 + c^2 * d );

      h(1,2) = ( 2.0 * a * b + a^2 * e ) ...
        * ( -2.0 * c * d + c^2 * ( r + s ) ) ...
        + ( 1.0 + a^2 * b ) ...
        * ( -12.0 * d + 4.0 * c * s -6.0 * c * r - 36.0 * c^2 ) ...
        + ( 2.0 * b + 4.0 * a * e + 6.0 * a^2 ) ...
        * ( 30.0 + c^2 * d );

      h(2,1) = h(1,2);

      h(2,2) = 2.0 * ( 2.0 * a * b + a^2 * e ) ...
        * ( -6.0 * c * d + c^2 * s ) ...
        + ( 1.0 + a^2 * b ) ...
        * ( 18.0 * d - 6.0 * c * s - 6.0 * c * s + 54.0 * c^2 ) ...
        + ( 2.0 * b + 2.0 * a * e + 2.0 * a * e + 6.0 * a^2 ) ...
        * ( 30.0 + c^2 * d );
    end
  end
end
