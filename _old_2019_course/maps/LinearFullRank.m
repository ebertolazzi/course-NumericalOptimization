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
        N = int32(2);
        M = int32(2);
      elseif nargin == 1
        N = varargin{1};
        M = varargin{1};
      elseif nargin > 2
        error('LinearFullRank: too many arguments in constructor');
      else
        N = varargin{1};
        M = varargin{2}; 
      end
      if ~isinteger(N) || ~isinteger(M)
        error('LinearFullRank: arguments must be integers, found %s and %s',class(N), class(M));
      end
      if N <= 1
        error('LinearFullRank: argument n must be an integer > 1, found %d',N);
      end
      if M < N
        error('LinearFullRank: argument m must be an integer >= n, found %d',M);
      end
      self@FunctionMap(int32(M),int32(N));
      self.exact_solutions = -ones(1,M); % f=m-n;
      self.guesses         = ones(1,N);
      
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate the entries (not squared).
      self.check_x(x);
      F    = zeros(self.M,1);
      for i = 1:self.M
          F(i) = -1;
          if i <= self.N
                F(i) = F(i) + x(i);
                for j=1:self.N
                    F(i) = F(i) -2/self.M*x(j);
                end   
          else
                for j=1:self.N
                    F(i) = F(i) -2/self.M*x(j);
                end
          end    
      end    
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      J    = zeros(self.M,self.N);
      for i = 1:self.M
          for j=i:self.N
              J(i,j) = -2/self.M;
              if j == i
                  J(i,j) = 1-2/self.M;
              end    
          end
     end   
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      T    = zeros(self.M,self.N,self.N);
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
