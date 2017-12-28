classdef Brown_bsf < FunctionND
  methods
    function self = Brown_bsf()
      arity = 2;
      self@FunctionND(int32(arity)) ;
    end
    
    function f = eval(self,x)
      % Evaluate Brown badly scaled 2D function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      X = squeeze(x(1,:,:)) ;
      Y = squeeze(x(2,:,:)) ;
      f1 = X - 10^6;
      f2 = Y -2*10^(-6);
      f3 = X.*Y - 2;
      f = f1^2 + f2^2 + f3^2;
    end
  end
end