classdef Easy_function_3D < FunctionND
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Easy_function_3D()
      arity = 2;
      self@FunctionND(int32(arity));
      self.exact_solutions = [ 0;.0 ];  % one known solution
      self.guesses         = [ 2; 3 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      % evaluate a simple (2D) function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      X = squeeze(x(1,:,:));
      Y = squeeze(x(2,:,:));
      f = (X.^2)+ abs(Y);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Use finite difference for grad and hessian
    function g = grad( self, x )
      g = self.FD_grad( x );
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function H = hessian( self, x )
      H = self.FD_hessian( x );
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
