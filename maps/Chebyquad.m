classdef Chebyquad < FunctionMap
  %
  %   Chebyquad function 
  %
  %   Fletcher, R. "Function minimization without evaluating derivatives -
  %   A review.", Comput. J. 8, 1965, 33-41.
  %
  %   see in reference test N.35
  %
  %   @article{More:1981,
  %     author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},
  %     title   = {Testing Unconstrained Optimization Software},
  %     journal = {ACM Trans. Math. Softw.},
  %     volume  = {7},
  %     number  = {1},
  %     year    = {1981},
  %     pages   = {17--41},
  %     doi     = {10.1145/355934.355936}
  %    }
  %
  %   Author: Eleonora Isotta - University of Trento
  
  
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Chebyquad( varargin )
      if nargin == 0
        N = 2;
        M = 2;
      elseif nargin == 1
        N = varargin{1};
        M = varargin{1};
      elseif nargin > 2
        error('Chebyquad: too many arguments in constructor');
      else
        N = varargin{1};
        M = varargin{2}; 
      end
      if N <= 1
        error('Chebyquad: argument n must be an integer > 1, found %d',N);
      end
      if M < N
        error('Chebyquad: argument m must be an integer >= n, found %d',M);
      end
      self@FunctionMap(int32(M),int32(N));
      self.exact_solutions = []; % no exact solutions;
                                 % f=0 for m=n, 1<=n<=7, and n=9
                                 % f=3.51687...10^(-3) for m=n=8,
                                 % f=6.50395...10^(-3) for m=n=10
      self.guesses         = 1/(N+1)*(1:1:N);
      
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate the entries (not squared).
      self.check_x(x);
      F    = zeros(self.M,1);
      t    = zeros(self.M,self.N);
      t(1,:) = ones(1,self.N);
      t(2,:) = x;
      F(1) = 1;
      F(2) = 1/double(self.N)*sum(t(2,:))+1/3;
      for i = 3:self.M
          t(i,:) = 2*x.*t(i-1,:)-t(i-2,:);
          F(i) = 1/double(self.N)*sum(t(i,:));
          if isinteger(i/2)
              F(i) = F(i) + 1/(i^2-1);
          end    
      end    
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      J    = zeros(self.M,self.N);
      t    = zeros(self.M,self.N);
      t(1,:) = ones(1,self.N);
      t(2,:) = x;
      J(2,:) = 1/double(self.N)*ones(1,self.N);
      for i = 3:self.M
          t(i,:) = 2*x.*t(i-1,:)-t(i-2,:);
          J(i,:) = 1/double(self.N)*(2*t(i-1,:) + 2*double(self.N)*x.*J(i-1,:) - double(self.N)*J(i-2,:));
      end   
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function t = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      T    = zeros(self.M,self.N,self.N);
      J    = zeros(self.M,self.N);
      t    = zeros(self.M,self.N);
      t(1,:) = ones(1,self.N);
      t(2,:) = x;
      J(2,:) = 1/double(self.N)*ones(1,self.N);
      for i = 3:self.M
          t(i,:) = 2*x.*t(i-1,:)-t(i-2,:);
          J(i,:) = 1/double(self.N)*(2*t(i-1,:) + 2*double(self.N)*x.*J(i-1,:) - double(self.N)*J(i-2,:));
          for j = 1 : self.N
              T(i,j,j) = 4*J(i-1,j) + 2*x(j)*T(i-1,j,j) - T(i-2,j,j);  
          end
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
