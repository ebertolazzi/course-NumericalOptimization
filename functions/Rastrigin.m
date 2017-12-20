classdef Rastrigin < FunctionND
  methods
    function self = Rastrigin()
      arity = 2;
      self@FunctionND(int32(arity)) ;
    end
    
    function f = eval(self,x)
      % evaluate Rosenbrock (2D) function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      X = squeeze(x(1,:,:)) ;
      Y = squeeze(x(2,:,:)) ;
      
      f = 20 + X.^2 - 10*cos(2*pi*X.^2) + Y.^2 - 10*cos(2*pi*Y.^2);
      
    end
  end
end