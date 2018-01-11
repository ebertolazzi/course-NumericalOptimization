classdef Himmelblau < FunctionMap
  %
  % The Himmelblau function
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
  %   David Himmelblau,
  %   Applied Nonlinear Programming,
  %   McGraw Hill, 1972,
  %   ISBN13: 978-0070289215,
  %   LC: T57.8.H55.
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

    function self = Himmelblau()
      self@FunctionMap(int32(2),int32(2)) ;
      self.exact_solutions = [ 3.0; 2.0];     % one known solution 
      self.guesses         = [ -1.3; 2.7];
    end

    function F = evalMap(self,x)
      % evaluate function
      self.check_x(x);
      X = x(1) ;
      Y = x(2) ;
      F = [ X^2 + Y - 11 ; X + Y^2 - 7 ] ;
    end

    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1) ;
      Y = x(2) ;
      J = [ 2*X, 1 ; 1, 2*Y ] ;
    end

    function T = tensor( self, x )
      % use analitic tensor of second derivative
      T = zeros(2,2,2);
      T(1,:,:) = [ 2, 0 ; 0, 0 ] ;
      T(2,:,:) = [ 0, 0 ; 0, 2 ] ;
    end

  end
end
