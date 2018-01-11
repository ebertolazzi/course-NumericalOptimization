classdef Leon < FunctionMap
  %
  % The Leon cubic valley function.
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %   A Leon,
  %   A Comparison of Eight Known Optimizing Procedures,
  %   in Recent Advances in Optimization Techniques,
  %   edited by Abraham Lavi, Thomas Vogl,
  %   Wiley, 1966.
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

    function self = Leon()
      self@FunctionMap(int32(2),int32(2)) ;
      self.exact_solutions = [ 1.0; 1.0  ];     % one known solution 
      self.guesses         = [ -1.2; -1.0 ];
    end

    function F = evalMap(self,x)
      % evaluate function
      self.check_x(x);
      X  = x(1) ;
      Y  = x(2) ;
      f1 = Y - X^3 ;
      f2 = 1.0 - X ;
      F  = [ sqrt(100.0) * f1 ; f2 ];
    end

    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1) ;
      J = [ -sqrt(100.0) * 3*X^2, sqrt(100.0) ; 0, -1 ] ; 
    end

    function T = tensor( self, x )
      % use analitic tensor of second derivative
      X        = x(1) ;
      T        = zeros(2,2,2);
      T(1,:,:) = [ -6*X, 0 ; 0, 0 ] ;
      T(2,:,:) = [ 0, 0 ; 0, 0 ] ;
    end

  end
end
