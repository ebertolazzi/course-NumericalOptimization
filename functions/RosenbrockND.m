classdef Rosenbrock < FunctionND % Funz_strana deriva dalla classe "astratta" Function1D
  methods
    function self = Rosenbrock()
    end
    
    function f = eval(self,x)
      n = length(x); % Domain space dimension;
      f = 0;
      for i = 1:n-1
        f = f + 100*(x(i+1)- x(i)^2)^2 + (1-x(i))^2;
      end
      
    end
  end
end