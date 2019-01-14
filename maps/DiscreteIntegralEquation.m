classdef DiscreteIntegralEquation < FunctionMap
  %
  %   Discrete integral equation function 
  %
  %   Morè, J. J., and Cosnard, M. Y. "Numerical solution of nonlinear
  %   equations". ACM Trans. Math. Softw 5, 1 (March 1979), 64-85.
  % 
  %
  %   see also in reference test N.29
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
  properties
      h % constant n-dependent
      t % vector h-dependent
  end
  %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = DiscreteIntegralEquation( varargin )
      if nargin == 0
        n = int32(2);
      elseif nargin == 1
        n = varargin{1};
      else
        error('DiscreteIntegralEquation: too many arguments in constructor');
      end
      if ~isinteger(n)
        error('DiscreteIntegralEquation: argument must be an integer, found %s',class(n));
      end
      if n <= 1
        error('DiscreteIntegralEquation: argument must be an integer > 1, found %d',n);
      end
      m = n;
      self@FunctionMap(int32(m),int32(n));
      self.h = 1/(n+1);
      self.t = h*(1:1:n);
      self.exact_solutions = zeros(m,0); % unknown solutions, f=0;
      self.guesses         = self.t.*(self.t-1);
      
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate the entries (not squared).
      self.check_x(x);
      F    = zeros(m,1);
      for i = 1:m
          F(i) = x(i);
          for j=1:n
              if j <= i
                    F(i) = F(i)+self.h/2*(1-self.t(i))*self.t(j)*(x(j)+self.t(j)+1)^3; 
              else
                    F(i) = F(i)+self.h/2*self.t(i)*(1-self.t(j))*(x(j)+self.t(j)+1)^3;
              end    
          end    
      end    
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      J    = zeros(m,n);
      for i = 1:m
          for j=i:n
              if j < i
                  J(i,j) = self.h/2*(1-self.t(i))*self.t(j)*3*(x(j)+self.t(j)+1)^2;
              elseif j == i
                  J(i,j) = 1+self.h/2*(1-self.t(i))*self.t(j)*3*(x(j)+self.t(j)+1)^2;
              else
                  J(i,j) = self.h/2*self.t(i)*(1-self.t(j))*3*(x(j)+self.t(j)+1)^2;
              end    
          end
     end   
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      T    = zeros(m,n,n);
      for i = 1:m
          for j = 1:n
              if j <= i
                  T(i,j,j) = self.h*(1-self.t(i))*self.t(j)*3*(x(j)+self.t(j)+1);
              else
                  T(i,j,j) = self.h*self.t(i)*(1-self.t(j))*3*(x(j)+self.t(j)+1);
              end
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
