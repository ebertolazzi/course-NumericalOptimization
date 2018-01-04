classdef (Abstract) FunctionND < handle

  properties (SetAccess = protected, Hidden = true)
    N               % arity of the function
    exact_solutions % matrix N x dim with all the known solutions. dim can be 0 if no known solution arer available 
    guesses         % matrix N x dim with suggested inital guess used for testing.
  end

  methods (Abstract)
    y = eval( self, x )
    J = grad( self, x )
    H = hessian( self, x )
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
      self.guesses         = zeros(N,0); % no inital guess
      self.exact_solutions = zeros(N,0); % no exact solutions
    end

    function n = num_guess( self )
      n = size(self.guesses,2);
    end

    function g = guess( self, idx )
      if ~isscalar(idx)
        error('FunctionND:guess, argument must be a scalar, found %s',class(idx));
      end
      if ~isinteger(idx)
        error('FunctionND:guess, argument must be an integer, found %s',class(idx));
      end
      if idx < 1 || idx > size(self.guesses,2)
        error('FunctionND:guess, argument must be an integer in [1,%d] found %s',size(self.guesses,2),idx);
      end
      g = self.guesses(:,idx);
    end

    function n = num_exact( self )
      n = size(self.exact_solutions,2);
    end

    function e = exact( self, idx )
      if ~isscalar(idx)
        error('FunctionND:exact, argument must be a scalar, found %s',class(idx));
      end
      if ~isinteger(idx)
        error('FunctionND:exact, argument must be an integer, found %s',class(idx));
      end
      if idx < 1 || idx > size(self.exact_solutions,2)
        error('FunctionND:exact, argument must be an integer in [1,%d] found %s',size(self.exact_solutions,2),idx);
      end
      e = self.exact_solutions(:,idx);
    end

    function N = arity( self )
      % return the number of arguments of the function 
      N = self.N ;
    end

    function check_x( self, x )
      szdim = size(x);
      if szdim ~= 2
        error('FunctionND, dimension of x = %%d expected 2\n',szdim);
      end
      n = size(x,1);
      m = size(x,2);
      if ~( n == self.N && m == 1 )
        error('FunctionND, size(x) = %d x %d, expected %d x 1\n',n,m,self.N);
      end 
    end

    function g = FD_grad( self, x )
      % finite difference approximation of the gradient
      h  = max(1,abs(x))*eps^(1/3);
      xp = x ;
      xm = x ;
      g  = zeros(1,self.N);
      for k=1:self.N
        xp(k) = x(k)+h(k);
        xm(k) = x(k)-h(k);
        fp    = self.eval(xp) ;
        fm    = self.eval(xm) ;
        g(k)  = (fp-fm)./(2*h(k)) ;
        xp(k) = x(k);
        xm(k) = x(k);
      end
    end

    function H = FD_hessian( self, x )
      % finite difference approximation of the hessian
      % Baseed on a code by Brendan C. Wood
      % Copyright (c) 2011, Brendan C. Wood <b.wood@unb.ca>
      H = zeros(self.N,self.N);
      h = max(1,abs(x))*eps^(1/3); %0.1*ones(3,1); % ricetta di cucina, attenzione x e' un vettore
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
      end

      % Make H symmetric numerically, this is not mandatory but could help
      H = 0.5*(H+H.') ;
    end

    function contour( self, xmm, ymm, fz, nc )
      % plot 2D contour of the function
      if self.N ~= 2
        error('FunctionND:contour can be used only for 2D functions');
      end
      NX    = 400 ;
      NY    = 400 ;
      x     = linspace(xmm(1),xmm(2),NX);
      y     = linspace(ymm(1),ymm(2),NY);
      [X,Y] = meshgrid(x,y);
      %
      % X and Y are matrices, to evaluate in [X(i,j);Y(i,j)]
      %
      Z = zeros(NX,NY) ;
      for i=1:NX
        for j=1:NY
         XY = [ X(i,j) ; Y(i,j) ] ;
         Z(i,j) = feval(fz,self.eval(XY)) ;
        end
      end
      contour(X,Y,Z,nc);
    end
  end
end