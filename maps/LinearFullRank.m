classdef LinearFullRank < FunctionMap
  %
  %   Linear function - full rank 
  %
  %   Unpublished.
  %
  %   see in reference test N.32
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
    function self = LinearFullRank( varargin )
      if nargin == 0
        n = int32(2);
        m = int32(2);
      elseif nargin == 1
        n = varargin{1};
        m = varargin{1};
      elseif nargin > 2
        error('LinearFullRank: too many arguments in constructor');
      else
        n = varargin{1};
        m = varargin{2}; 
      end
      if ~isinteger(n) || ~isinteger(m)
        error('LinearFullRank: arguments must be integers, found %s and %s',class(n), class(m));
      end
      if n <= 1
        error('LinearFullRank: argument n must be an integer > 1, found %d',n);
      end
      if m < n
        error('LinearFullRank: argument m must be an integer >= n, found %d',m);
      end
      self@FunctionMap(int32(m),int32(n));
      self.exact_solutions = -ones(1,m); % f=m-n;
      self.guesses         = ones(1,n);
      
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate the entries (not squared).
      self.check_x(x);
      F    = zeros(m,1);
      for i = 1:m
          F(i) = -1;
          if i <= n
                F(i) = F(i) + x(i);
                for j=1:n
                    F(i) = F(i) -2/m*x(j);
                end   
          else
                for j=1:n
                    F(i) = F(i) -2/m*x(j);
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
              J(i,j) = -2/m;
              if j == i
                  J(i,j) = 1-2/m;
              end    
          end
     end   
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      T    = zeros(m,n,n);
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
