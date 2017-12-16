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
    use_FD % flag true = use finite difference in evaluation of first and second derivative
  end
  
  methods

    function obj = Function1Dcut( object_function_ND, x0, d, varargin )
      %% Constructor
      % Given the object object_functionND which represents
      % a N-dimensional function build the cut function $g(\alpha)$ 
      % in the d direction starting from $x_0$.

      obj.funND  = object_function_ND ;

      % check arguments
      obj.x0     = x0 ;
      obj.d      = d ;
      obj.use_FD = true ;
      if nargin >= 4 ; obj.use_FD = varargin{1} ; end
    end
    
    function res = eval( self, alpha )
      % This method evaluate the function in a point of the "line of cut".
      % $\alpha$ represents the coordinate along the line.
      res = self.funND.eval( self.x0 + alpha * self.d );
    end

    function res = ANALITIC_eval_D( self, alpha )
      % This method evaluate the first derivative of the cut function.
      % Used only for comparison purpose
      res = dot( self.funND.grad( self.x0 + alpha * self.d ), self.d ) ;
    end
    
    function res = ANALITIC_eval_DD( self, alpha )
      % This method evaluate the first derivative of the cut function.
      % Used only for comparison purpose
      res = dot( self.funND.hessian( self.x0 + alpha * self.d )*self.d', self.d ) ; 
    end
    
  end
end
