classdef Quadratic2D < FunctionND
  methods
    function self = Quadratic2D()
      self@FunctionND(int32(2)) ;
    end
    
    function f = eval(self,x)
      X = squeeze(x(1,:,:)) ;
      Y = squeeze(x(2,:,:)) ;
      f = 100*X.^2+Y.^2;
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