classdef Function1Dcut < Function1D
  % This class take a N_D function, a starting point and a search
  % direction, and "cuts" the original function in that direction.
  % In this way we "reduce" the function to a 1-D one.
  % First and second derivative (eval_D and evald_DD) by finite difference are
  % implemented in the base abstract class Function1D

  properties (SetAccess = private, Hidden = true)
    funND  % Original input function
    x0     % Initial point (n dimensional array, which express a point on the domain of the function)
    d      % Search direction (n dimensional array, which express the direction in which we have to search for the minimum)
    FD_D   % flag true = use finite difference in evaluation of first derivative
    FD_DD  % flag true = use finite difference in evaluation of second derivative
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Function1Dcut( object_function_ND, x0, d )
      %% Constructor
      % Given the object object_functionND which represents
      % a N-dimensional function build the cut function $g(\alpha)$
      % in the d direction starting from $x_0$.

      self.funND = object_function_ND;

      if norm(d,inf) == 0
        error('Function1Dcut, bad direction d == 0\n');
      end

      % check arguments
      self.x0    = x0;
      self.d     = d;
      self.FD_D  = true;
      self.FD_DD = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function use_FD_D( self )
      self.FD_D = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function no_FD_D( self )
      self.FD_D = false;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function use_FD_DD( self )
      self.FD_D = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function no_FD_DD( self )
      self.FD_D = false;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = eval( self, alpha )
      % This method evaluate the function in a point of the "line of cut".
      % $\alpha$ represents the coordinate along the line.
      res = self.funND.eval( self.x0 + alpha * self.d );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = eval_D( self, alpha )
      if self.FD_D
        res = self.FD_eval_D( alpha );
      else
        g   = self.funND.grad( self.x0 + alpha * self.d );
        res = dot( g.', self.d );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = eval_DD( self, alpha )
      if self.FD_DD
        res = self.FD_eval_DD( alpha );
      else
        H   = self.funND.hessian( self.x0 + alpha * self.d );
        res = dot( H * self.d, self.d );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g] = eval_FG( self, x )
      f = self.eval(x);
      g = self.eval_D(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [f,g,H] = eval_FGH( self, x )
      f = self.eval(x);
      g = self.eval_D(x);
      H = self.eval_DD(x);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
