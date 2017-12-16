classdef Easy_function_4D < FunctionND % Funz_strana deriva dalla classe "astratta" Function1D
  methods
    function self = Easy_function_4D()
    end
    
    function fxyz = eval(self,x)

      fxyz = x(1)^2 + x(2) + x(3);
      
    end
  end
end