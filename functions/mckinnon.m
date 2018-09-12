classdef mckinnon < FunctionND
  %
  % The McKinnon function
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
  %   Ken McKinnon,
  %   Convergence of the Nelder-Mead simplex method to a nonstationary point,
  %   SIAM Journal on Optimization,
  %   Volume 9, Number 1, 1998, pages 148-158.
  %
  %  Discussion:
  %
  %   This function has a global minimizer:
  %
  %     X* = ( 0.0, -1.0 ), F(X*) = 0.
  %
  %   There are three parameters, TAU, THETA and PHI.
  %
  %   1 < TAU, then F is strictly convex.
  %            and F has continuous first derivatives.
  %   2 < TAU, then F has continuous second derivatives.
  %   3 < TAU, then F has continuous third derivatives.
  %
  %   However, this function can cause the Nelder-Mead optimization
  %   algorithm to "converge" to a point which is not the minimizer
  %   of the function F.
  %
  %   Sample parameter values which cause problems for Nelder-Mead
  %   include:
  %
  %     TAU = 1, THETA = 15, PHI =  10;
  %     TAU = 2, THETA =  6, PHI =  60;
  %     TAU = 3, THETA =  6, PHI = 400;
  %
  %   To get the bad behavior, we also assume the initial simplex has the form
  %
  %     X1 = (0,0),
  %     X2 = (1,1),
  %     X3 = (A,B),
  %
  %   where
  %
  %     A = (1+sqrt(33))/8 =  0.84307...
  %     B = (1-sqrt(33))/8 = -0.59307...
  %
  % Author:
  %
  %   Enrico Bertolazzi
  %   Dipartimento di Ingegneria Indutriale
  %   Universita` degli Studi di Trento
  %   email: enrico.bertolazzi@unitn.it
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %

  properties (SetAccess = protected, Hidden = true)
    tau
    theta
    phi
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = mckinnon()
      self@FunctionND(int32(2));
      self.exact_solutions = [ 0; -1];     % one known solution
      self.guesses         = [ 1; 1];
      self.tau   = 2.0;
      self.theta = 6.0;
      self.phi   = 60.0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      self.check_x(x);
      if x(1) <= 0.0
        t2 = (-x(1)) ^ self.tau;
        f  = t2 * self.theta * self.phi + x(2) * (1 + x(2));
      else
        t1 = x(1) ^ self.tau;
        f  = self.theta * t1 + x(2) * (1 + x(2));
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 2 );
      if x(1) <= 0.0
        t2   = (-x(1)) ^ self.tau;
        g(1) = self.theta * self.phi * t2 * self.tau / x(1);
        g(2) = 1 + (2 * x(2));
      else
       t1   = x(1) ^ self.tau;
       g(1) = self.theta * t1 * self.tau / x(1);
       g(2) = 1 + (2 * x(2));
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( 2, 2 );
      if ( x(1) <= 0.0 )
        t1 = self.theta * self.phi;
        t2 = (-x(1)) ^ self.tau;
        t3 = self.tau ^ 2;
        t5 = x(1) ^ 2;
        h(1,1) = (t3 - tau) * t1 * t2/ t5;
        h(1,2) = 0;
        h(2,1) = 0;
        h(2,2) = 2;
      else
        t1 = x(1) ^ self.tau;
        t2 = self.theta * t1;
        t3 = self.tau ^ 2;
        t4 = x(1) ^ 2;
        h(1,1) = (t3 - tau) * t2./t4;
        h(1,2) = 0;
        h(2,1) = 0;
        h(2,2) = 2;
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
