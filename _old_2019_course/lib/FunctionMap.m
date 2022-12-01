classdef (Abstract) FunctionMap < FunctionND
% Function of the type "sum of squares"
%
% Map a sum of squares function from $f_i(x)$ to $\frac{1}{2}\sum_{i=1}^{m}f_i^2(x)$.


  properties (SetAccess = protected, Hidden = true)
    M % number of components of the Map
  end

  methods (Abstract)
    F = evalMap( self, x )
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = FunctionMap( N, M )
      % Constructor of base (abstract) class
      self@FunctionND(N); % call the contructor of the superclass
      if ~isinteger(M)
        error('FunctionND: argument must be an integer, found %s',class(M));
      end
      if M <= 0
        error('FunctionND: argument must be a positive integer, found %d',M);
      end
      self.M = M;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = is_a_map( self )
      ok = true;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function m = dimension( self )
      % return the number of components of the map
      m = self.M;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = FD_jacobian( self, x )
      % finite difference approximation of the jacobian
      self.check_x(x);
      h  = max(1,abs(x))*eps^(1/3);
      xp = x;
      xm = x;
      J  = zeros(self.M,self.N);
      for k=1:self.N
        xp(k)  = x(k)+h(k);
        xm(k)  = x(k)-h(k);
        fp     = self.evalMap(xp);
        fm     = self.evalMap(xm);
        J(:,k) = (fp-fm)./(2*h(k));
        xp(k)  = x(k);
        xm(k)  = x(k);
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = FD_tensor( self, x )
      % finite difference approximation of the second derivative
      % T(k,:,:) is the hessian of the k-th components of the map
      self.check_x(x);
      T = zeros(self.M,self.N,self.N);
      h = max(1,abs(x))*eps^(1/3); %0.1*ones(3,1); % ricetta di cucina, attenzione x e' un vettore
      % for each dimension of objective function
      for i=1:self.N
        % derivative at first point (left)
        x1    = x;
        x1(i) = x(i) - h(i);
        jf1   = self.jacobian(x1);

        % derivative at second point (right)
        x2    = x;
        x2(i) = x(i) + h(i);
        jf2   = self.jacobian(x2);

        % differentiate between the two derivatives
        j2f = (jf2-jf1) ./ (2*h(i)); % M x N

        % assign as column i of Hessian
        T(:,i,:) = j2f;
      end

      % Make T symmetric numerically, this is not mandatory but could help
      for k=1:self.M
        Ttemp    = squeeze(T(k,:,:));     % added this temporary variable to use the .' operator
        T(k,:,:) = 0.5*(Ttemp + Ttemp.');
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % rewrite of gradient and hessian of the sum of the squares of the map
    % Virtual methods definitions
    function f = eval( self, x )
      % compute the function as the sum of the squares of the components of the map
      self.check_x(x);
      F = self.evalMap( x );
      f = 0.5*dot(F,F);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function g = grad( self, x )
      % gradient of sum f_k(x)^2 --> sum f_k(x) grad_k(x) = 2*F(x)^T J(x)
      self.check_x(x);
      F = self.evalMap(x);
      J = self.jacobian(x);
      g = F.' * J;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function H = hessian( self, x )
      % hessian of sum f_k(x)^2 --> 2*(J^T(x) *J(x) + sum f_k(x) Hessian_k(x))
      self.check_x(x);
      F = self.evalMap( x );
      J = self.jacobian( x );
      T = self.tensor( x );
      H = J.' * J;
      for k=1:self.M
        H = H + F(k) * squeeze(T(k,:,:));
      end
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
