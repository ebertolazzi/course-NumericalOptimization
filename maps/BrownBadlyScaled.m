classdef BrownBadlyScaled < FunctionMap
  methods
    function self = BrownBadlyScaled()
      self@FunctionMap(int32(2),int32(3));
      self.exact_solutions = [ 1.0; 1.0 ];     % one known solution
      self.guesses         = [ 1e6; 2e6 ];
    end

    function F = evalMap(self,x)
      % Evaluate Brown badly scaled 2D function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      self.check_x(x);
      X = x(1);
      Y = x(2);
      F = [ X - 1e6; Y - 2e-6; X*Y - 2 ];
    end

    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1);
      Y = x(2);
      J = [ 1, 0; 0, 1; Y, X ];
    end

    function T = tensor( self, x )
      % use analitic tensor of second derivative
      T = zeros(3,2,2);
      T(1,:,:) = [ 0, 0; 0, 0 ];
      T(2,:,:) = [ 0, 0; 0, 0 ];
      T(3,:,:) = [ 0, 1; 1, 0 ];
    end

  end
end
