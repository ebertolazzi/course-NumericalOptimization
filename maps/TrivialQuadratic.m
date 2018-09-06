classdef TrivialQuadratic < FunctionMap
  %
  % Function de mierda
  %

  methods

    function self = TrivialQuadratic()
      self@FunctionMap(int32(2),int32(2));
      self.exact_solutions = [0;0];   % one known solution
      self.guesses         = [-1;1]; % one guess
    end

    function F = evalMap(self,x)
      % evaluate quadratic (2D) function.
      X = x(1);
      Y = x(2);
      F = [ X; Y ];
    end

    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1);
      J = [ 1 , 0; 0 , 1 ];
    end

    function T = tensor( self, x )
      % use analitic tensor of second derivative
      T        = zeros(2,2,2);
      T(1,:,:) = [ 0, 0; 0, 0 ];
      T(2,:,:) = [ 0, 0; 0, 0 ];
    end

  end
end
