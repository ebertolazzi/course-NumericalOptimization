classdef Rosenbrock < FunctionND % Funz_strana deriva dalla classe "astratta" Function1D
  methods
    function self = Rosenbrock()
      self@FunctionND(int32(2)) ;
    end
    
    function f = eval(self,x)
      % evaluate Rosenbrock (2D) function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      X = squeeze(x(1,:,:)) ;
      Y = squeeze(x(2,:,:)) ;
      f = 100*(Y-X.^2).^2 + (1-X).^2;
    end
  end
end