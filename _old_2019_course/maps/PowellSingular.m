classdef PowellSingular < FunctionMap
  %
  % The Powerll Singular Function
  %
  % Reference:
  %
  %   Powell, M. J. D.
  %   An Iterative Method for Finding Stationary Values of a Function of Several Variables
  %   The Computer Journal, vol 5, n.2, pp. 147-151, 1962,
  %
  % Richard Brent,
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
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = PowellSingular()
      self@FunctionMap(int32(4),int32(4));
      self.exact_solutions = [ 101; 10; 0; 0 ]; % one known solution
      self.guesses         = [ 0.0; 0.0; 0.0; 0.0 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate function
      self.check_x(x);
      F = [ x(1) + 10*x(2); ...
            sqrt(5)*(x(3)-x(4)); ...
            (x(2)-2*x(3))^2; ...
            sqrt(10)*(x(1)-x(4))^2 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      t1 = x(2)-2*x(3);
      t2 = sqrt(10)*(x(1)-x(4));
      J = [ 1,    10,   0,       0; ...
            0,    0,    sqrt(5), -sqrt(5); ...
            0,    2*t1, -4*t1,   0; ...
            2*t2, 0,    0,       -2*t2 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analitic tensor of second derivative
      T        = zeros(4,4,4);
      T(3,:,:) = [ 0,  0,  0,  0; ...
                   0,  2, -4,  0; ...
                   0, -4,  8,  0; ...
                   0,  0,  0,  0 ];
      T(4,:,:) = 2*sqrt(10)*[ 1,  0,  0, -1; ...
                              0,  0,  0,  0; ...
                              0,  0,  0,  0; ...
                             -1,  0,  0,  1 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
