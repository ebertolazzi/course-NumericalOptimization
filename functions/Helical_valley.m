classdef Helical_valley < FunctionND
  methods
    function self = Helical_valley()
      arity = 3;
      self@FunctionND(int32(arity)) ;
    end
    
    function f = eval(self,x)
      % Evaluate Brown badly scaled 2D function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      X1 = squeeze(x(1,:,:)) ;
      X2 = squeeze(x(2,:,:)) ;
      X3 = squeeze(x(3,:,:)) ;
      
      f_theta  = @(X1,X2) 1/(2*pi)*piecewise(X1 > 0, atan(X2/X1), X1 < 0  , atan(X2/X1)+0.5   ); 
      
      f1 = 10*( X3    - 10*f_theta(X1,X2) );
      f2 = 10*( sqrt( X1.^2 + X2.^2)  - 1 );
      f3 = X3;
      f = f1^2 + f2^2 + f3^2;
      
    end
  end
end