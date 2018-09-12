classdef Rastrigin < FunctionND
  methods
    function self = Rastrigin()
      arity = 2;
      self@FunctionND(int32(arity));
      self.exact_solutions = [ 0;0];
      self.guesses         = [ 2;3];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function f = eval(self,x)
      % evaluate Rastrigin (2D) function.
      % if x is a 2 by m matrix return m values in a row vector.
      % if x is a 2 by m x n matrix return m x n values in a matrix vector.
      X = squeeze(x(1,:,:));
      Y = squeeze(x(2,:,:));
      f = 20 + X.^2 + Y.^2 - 10* (cos(2*pi*X.^2)+cos(2*pi*Y.^2));
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad(self,x)
      % evaluate Rastrigin (2D) function.
      g(1) = 2*x(1)+40*pi*x(1)*sin(2*pi*x(1)^2);
      g(2) = 2*x(2)+40*pi*x(2)*sin(2*pi*x(2)^2);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Use finite difference for hessian
    function h = hessian( self, x )
      h = self.FD_hessian( self, x );
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
