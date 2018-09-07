classdef Box3 < FunctionMap
  %
  % The Box 3-dimensional function.
  %
  % Reference:
  %   Box, M. J.
  %   A Comparison of Several Current Optimization Methods, and the use of Transformations in Constrained Problems
  %   The Computer Journal, vol 9, n.1, pp. 67-77, 1966,
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
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Box3()
      self@FunctionMap(int32(3),int32(10));
      self.exact_solutions = [ 1.0; 10.0; 1.0];     % one known solution
      self.guesses         = [ 0.0; 10.0; 5.0 ];
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate function
      self.check_x(x);
      F = zeros(10,1);
      for i = 1:10
        c    = - i / 10.0;
        F(i) = exp( c * x(1) ) - exp( c * x(2) ) - x(3) * ( exp( c ) - exp( 10.0 * c ) );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      J = zeros(10,3);
      for i = 1:10
        c      = - i / 10.0;
        J(i,:) = [ c * exp( c * x(1) ), - c * exp( c * x(2) ),  exp( 10.0 * c )  - exp( c ) ];
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analitic tensor of second derivative
      T = zeros(10,3,3);
      for i = 1:10
        c        = - i / 10.0;
        T(i,:,:) = c^2 * [ exp( c * x(1) ), 0, 0;  ...
                           0, -exp( c * x(2) ), 0; ...
                           0, 0, 0 ];
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
