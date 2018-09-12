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
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Box3()
      self@FunctionMap(int32(3),int32(10));
      self.exact_solutions = [ 1.0; 10.0; 1.0];     % one known solution
      self.guesses         = [ 0.0; 10.0; 5.0 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap( self, xx )
      % evaluate function
      self.check_x( xx );
      x1 = xx(1);
      x2 = xx(2);
      x3 = xx(3);
      F  = zeros(10,1);
      for i = 1:10
        c    = -i/10;
        F(i) = exp(c*x1) - exp(c*x2) - x3*( exp(c) - exp(10*c) );
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, xx )
      % use analitic jacobian
      self.check_x( xx );
      x1 = xx(1);
      x2 = xx(2);
      x3 = xx(3);
      J  = zeros(10,3);
      for i = 1:10
        c      = -i/10;
        J(i,:) = [ c * exp(c*x1), -c * exp(c*x2), exp(10*c) - exp(c) ];
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, xx )
      self.check_x( xx );
      x1 = xx(1);
      x2 = xx(2);
      x3 = xx(3);
      % use analitic tensor of second derivative
      T = zeros(10,3,3);
      for i = 1:10
        c        = - i / 10.0;
        T(i,:,:) = c^2 * [ exp(c*x1),          0, 0; ...
                           0,         -exp(c*x2), 0; ...
                           0,                  0, 0 ];
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
