classdef Quadratic2D < FunctionND
  methods
    function self = Quadratic2D()
      self@FunctionND(int32(2));
      self.exact_solutions = [ 0;0];
      self.guesses         = [ 2;3];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      X = squeeze(x(1,:,:));
      Y = squeeze(x(2,:,:));
      f = 100*X.^2+Y.^2;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Use finite difference for grad and hessian
    function g = grad( self, x )
      g = self.FD_grad( x );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function h = hessian( self, x )
      h = self.FD_hessian( x );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = eval_FG( self, x )
      f = self.eval(x);
      g = self.grad(x);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g,H] = eval_FGH( self, x )
      f = self.eval(x);
      g = self.grad(x);
      H = self.hessian(x);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
