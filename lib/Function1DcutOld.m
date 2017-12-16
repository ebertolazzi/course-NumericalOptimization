classdef Function1Dcut < handle
  % This class take a N_D function, a starting point and a search
  % direction, and "cuts" the original function in that direction.
  % In this way we "reduce" the function to a 1-D one. 
  % Now we can easilly evaluate the first and secon derivatives of the
  % reduced function.
  
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
      res = dot( self.funND.grad( self.x0 + alpha * self.d ), self.d ) ;
    end
    
    function res = ANALITIC_eval_DD( self, alpha )
      % This method evaluate the first derivative of the cut function.
      res = dot( self.funND.hessian( self.x0 + alpha * self.d )*self.d', self.d ) ; % valuto l'hessiano nel punto nuvo più la direzione
    end
    
    function res = FD_eval_D( self, alpha )
      % This method evaluate the first derivative of the cut function by finite difference.
      h   = (1+abs(alpha))*eps^(1/3);
      fp  = self.funND.eval( self.x0 + (alpha+h) * self.d);
      fm  = self.funND.eval( self.x0 + (alpha-h) * self.d);
      res = (fp -fm)/(2*h);
    end

    function res = FD_eval_DD( self, alpha )
      % This method evaluate the second derivative of the cut function by finite difference.
      h   = (1+abs(alpha))*eps^(1/3);
      fp  = self.eval(alpha+h) ;
      fm  = self.eval(alpha-h) ;
      fc  = self.eval(alpha) ;
      res = (fp+fm-2*fc)./(h.^2) ;
    end
  end
end
