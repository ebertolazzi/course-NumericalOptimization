classdef Schwefel < FunctionND % Funz_strana deriva dalla classe "astratta" Function1D
  methods
    function self = Schwefel()
    end
    
    function f = eval(self,x)
      n = length(x); % Domain space dimension;
      f = 0;
      for i = 1:n
        f = f - x(i)*sin( sqrt(abs(x(i))) );
      end
    end
  end
end