classdef Brown_bsf < FunctionND
  %
  % The Brown Badly Scaled Function (unpublished).
  %
  % Reference:
  %	  This function is number 4 in the paper:
  %
  %   Morè J., Garbow S. and Hillstrom E. 
  %   Testing Uncostrained Optimization Software,
  % 
  %   , 1981,
  %   ISBN: 
  %   LC:
  %
  % Author: Davide Vignotto
  %   

  methods
    function self = Brown_bsf()
	  self@FunctionND(int32(2)) ;
      self.exact_solutions = [ 10^6 ; 2*10^(-6) ]; % one known solution 
      self.guesses         = [ 1; 1] ;
    end
    
    function f = eval(self,x)
      % Evaluate Brown badly scaled 2D function.
	    self.check_x(x);
      f1 = x(1) - 10^6;
      f2 = x(2) -2*10^(-6);
      f3 = x(1)*x(2) - 2;
      f = f1^2 + f2^2 + f3^2;
    end
    

    function g = grad( self, x )
      %g = self.FD_grad(x);
      %return
	    % use analitic gradient
      self.check_x(x);
      g  = zeros ( 1, 2 );
	  f1 = x(1) - 10^6;
      f2 = x(2) -2*10^(-6);
      f3 = x(1)*x(2) - 2;
	  g(1) = 2*f1 + 2*f3*x(2);
	  g(2) = 2*f2 + 2*f3*x(1);
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h  = zeros ( 2, 2 );
      h(1,1) = 2*(1 + x(2)^2);
      h(1,2) = 4*(x(1)*x(2) - 1);
      h(2,1) = h(1,2);
      h(2,2) =2*(1 + x(1)^2);
    end
 
  end
end