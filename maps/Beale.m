classdef Beale < FunctionMap
  %
  % The Beale function.
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
  %   Evelyn Beale,
  %   On an Iterative Method for Finding a Local Minimum of a Function of More than One Variable,
  %   Technical Report 25,
  %   Statistical Techniques Research Group,
  %   Princeton University, 1958.
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

  properties (SetAccess = private, Hidden = true)
    c
  end

  methods

    function self = Beale()
      self@FunctionMap(int32(2),int32(3)) ;
      self.exact_solutions = [ 1.0; 1.0  ];     % one known solution 
      self.guesses         = [ -1.2; -1.0 ];
      self.c = [ 1.5, 2.25, 2.625 ] ;
    end

    function F = evalMap(self,x)
      % evaluate function
      self.check_x(x);
      X = x(1) ;
      Y = x(2) ;
      F = [ self.c(1) - X*(1-Y) ; self.c(2) - X*(1-Y^2) ; self.c(3) - X*(1-Y^3) ] ;
    end

    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1) ;
      Y = x(2) ;
      J = [ Y-1, X ; Y^2-1, 2*X*Y ; Y^3-1, 3*X*Y^2 ] ; 
    end

    function T = tensor( self, x )
      % use analitic tensor of second derivative
      X = x(1) ;
      Y = x(2) ;
      T        = zeros(2,2,2);
      T(1,:,:) = [ 0, 1 ; 1, 0 ] ;
      T(2,:,:) = [ 0, 2*Y ; 2*Y, 2*X ] ;
      T(3,:,:) = [ 0, 3*Y^2 ; 3*Y^2, 6*X*Y ] ;
    end

  end
end
