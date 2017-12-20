classdef Easy_function_3D < FunctionND
  methods
    function self = Easy_function_3D()
      arity = 2;
      self@FunctionND(int32(arity)) ;
    end
    
    function f = eval(self,x)
      % evaluate a simple (2D) function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      X = squeeze(x(1,:,:)) ;
      Y = squeeze(x(2,:,:)) ;
      f = (X.^2)+ abs(Y);
    end
  end
end