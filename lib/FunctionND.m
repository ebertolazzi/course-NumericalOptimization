classdef (Abstract) FunctionND < handle

  properties (SetAccess = private, Hidden = true)
    N  % arity of the function
  end

  methods (Abstract)
    y = eval( self, x )
  end
  
  methods

    function self = FunctionND( N )
      % Constructor of base (abstract) class
      if ~isinteger(N)
        error('FunctionND: argument must be an integer, found %s',class(N));
      end
      if N <= 0
        error('FunctionND: argument must be a positive integer, found %d',N);
      end
      self.N = N ;
    end

    function N = arity( self )
      N = self.N ;
    end

    function g = grad( self, x )
      % finite difference approximation of the gradient
      h  = (1+abs(x))*eps^(1/3); % ricetta di cucina, attenzione x e' un vettore
      xp = x ;
      xm = x ;
      for k=1:self.N
        xp(k) = x(k)+h(k);
        xm(k) = x(k)-h(k);
        fp    = self.eval(xp) ;
        fm    = self.eval(xm) ;
        g(k)  = (fp-fm)./(2*h(k)) ; %#ok<AGROW> % vectorial expression g(i) = (fp(i)-fm(i))/(2*h(i))
        xp(k) = x(k);
        xm(k) = x(k);
      end
    end
    
    function H = hessian( self, x )
      % finite difference approximation of the hessian
      % Baseed on a code by Brendan C. Wood
      % Copyright (c) 2011, Brendan C. Wood <b.wood@unb.ca>
      H = zeros(self.N,self.N);
      h  = (1+abs(x))*eps^(1/3); %0.1*ones(3,1); % ricetta di cucina, attenzione x e' un vettore
      % for each dimension of objective function
      for i=1:self.N
        % derivative at first point (left)
        x1    = x ;
        x1(i) = x(i) - h(i) ;
        df1   = self.grad(x1);
        
        % derivative at second point (right)
        x2    = x;
        x2(i) = x(i) + h(i);
        df2   = self.grad(x2);
        
        % differentiate between the two derivatives
        d2f = (df2-df1) ./ (2*h(i));
        
        % assign as column i of Hessian
        H(:,i) = d2f;
        
        % Make H symmetric numerically, this is not mandatory but could
        % help
        H = 0.5*(H+H.') ;
      end
    end
  end
end