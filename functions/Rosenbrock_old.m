classdef Rosenbrock_old < FunctionND
  methods
    function self = Rosenbrock_old()
      arity = 2;
      self@FunctionND(int32(arity)) ;
    end
    
    function f = eval(self,x)
      % evaluate Rosenbrock (2D) function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      X = squeeze(x(1,:,:)) ;
      Y = squeeze(x(2,:,:)) ;
      f = 100*(Y-X.^2).^2 + (1-X).^2;
    end
    
    % Use finite difference for grad and hessian
    function g = grad( self, x )
      g = self.FD_grad( self, x );
    end

    function h = hessian( self, x )
      h = self.FD_hessian( self, x );
    end
  end
end