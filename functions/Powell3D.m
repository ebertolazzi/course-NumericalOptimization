classdef Powell3D < FunctionND
  %
  % The Bohachevsky Function #1. (function for problem 38)
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
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Powell3D()
      self@FunctionND(int32(3));
      self.exact_solutions = ones(3,1); % one known solution
      self.guesses         = [ 0.0; 1.0; 2.0 ];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,X)
      % evaluate function
      self.check_x(X);

      x = X(1);
      y = X(2);
      z = X(3);

      if ( y == 0.0 )
        term = 0.0;
      else
        %arg  = ( x + 2.0 * y + z ) / y;
        arg  = 2.0 + ( x + z ) / y;
        term = exp( - arg^2 );
      end

      f = 3.0 ...
        - 1.0 / ( 1.0 + ( x - y )^2 ) ...
        - sin ( 0.5 * pi * y * z ) ...
        - term;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, X )
      % use analitic gradient
      self.check_x(X);
      g = zeros( 1, 3 );
      x = X(1);
      y = X(2);
      z = X(3);

      tmp1 = 2*(x-y) / (1+(x-y)^2)^2;
      tmp2 = (pi/2)*cos( (pi/2)*y*z );

      g(1) = tmp1;
      g(2) = -tmp1 - tmp2 * z;
      g(3) =       - tmp2 * y;

      if ( y ~= 0.0 )
        tmp1 = (2*y+x+z)/y;
        tmp2 = 2*(2*y+x+z)*exp(-tmp1^2)/y^2;
        g(1) = g(1) - tmp2;
        g(2) = g(2) + tmp2*(x+z)/y;
        g(3) = g(3) - tmp2;
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, X )
      % use analitic hessian
      self.check_x(X);
      h = zeros ( 3, 3 );
      x = X(1);
      y = X(2);
      z = X(3);
      xy2  = (x-y)^2;
      tmp1 = 2-6*xy2;
      tmp2 = (1+xy2)^3;
      S    = sin((pi/2)*y*z);
      C    = cos((pi/2)*y*z);

      h(1,1) = tmp1/tmp2;
      h(1,2) = -h(1,1);
      h(1,3) = 0;

      h(2,2) = z^2*(pi/2)^2*S+(2-6*xy2)/tmp2;
      h(2,3) = (pi/4)*(pi*y*z*S-2*C);

      h(3,3) = (pi*y/2)^2*S;

      if ( y ~= 0.0 )
        t1  = (x + z);
        t2  = t1 + 2 * y;
        t3  = 1 / y;
        t4  = t3 ^ 2;
        t5  = t4 ^ 2;
        t6  = y ^ 2;
        t7  = 4 * y;
        t8  = t1 ^ 2;
        t9  = 3 * z;
        t2  = (0.4e1 * exp(-(t2 ^ 2 * t4)));
        t10 = t2 * ((t1 * (t7 + t1)) + 0.7e1 / 0.2e1 * t6) * t5;
        t3  = t2 * (t1 * t8 - y * t6 + t1 * y * (4 * t1 + 3 * y)) * t3 * t5;
        tmp = [ t10, ...
                -t3, ...
                t10, ...
                t2 * (y * ((0.5e1 / 0.2e1 * y + (4 * z)) * z - (2 * t6)) + (z ^ 3) + (((t7 + t9 + x) * x) + ((t9 + 8 * y) * z) + 0.5e1 / 0.2e1 * t6) * x) * t1 * t4 * t5, ...
                -t3, ...
                t10 ];
        h(1,1) = h(1,1) + tmp(1);
        h(1,2) = h(1,2) + tmp(2);
        h(1,3) = h(1,3) + tmp(3);
        h(2,2) = h(2,2) + tmp(4);
        h(2,3) = h(2,3) + tmp(5);
        h(3,3) = h(3,3) + tmp(6);
      end
      h(2,1) = h(1,2);
      h(3,1) = h(1,3);
      h(3,2) = h(2,3);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = eval_FG( self, x )
      f = self.eval(x);
      g = self.grad(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g,H] = eval_FGH( self, x )
      f = self.eval(x);
      g = self.grad(x);
      H = self.hessian(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
