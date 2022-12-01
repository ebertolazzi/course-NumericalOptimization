classdef Trigonometric < FunctionMap
  %
  %   Trigonomertic function 
  %
  %   Spedicato, E. "Computational experience with quasi-Newton algorithms
  %   for minimization problems of moderately large size". Rep. CISE-N-175, 
  %   Segrate (Milano), 1975.
  % 
  %
  %   see also in reference test N.26
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
    function self = Trigonometric( varargin )
      if nargin == 0
        N = int32(2);
      elseif nargin == 1
        N = varargin{1};
      else
        error('Trigonometric: too many arguments in constructor');
      end
      if ~isinteger(N)
        error('Trigonometric: argument must be an integer, found %s',class(N));
      end
      if N <= 1
        error('Trigonometric: argument must be an integer > 1, found %d',N);
      end
      M = N;
      self@FunctionMap(int32(M),int32(N)); 
      self.exact_solutions = zeros(M,0); % unknown solutions, f=0;
      self.guesses         = 1/double(N)*ones(1,N);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate the entries (not squared).
      self.check_x(x);
      F    = zeros(self.M,1);
      for i = 1:self.M
          F(i) = self.N + i*(1-cos(x(i))) - sin(x(i));
          for j=1:self.N
             F(i) = F(i) - cos(x(j)); 
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
              J(i,j) = sin(x(j));
              if j == i
                  J(i,j) = J(i,j) + i*sin(x(i)) - cos(x(i));
              end    
          end
     end   
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      T    = zeros(self.M,self.N,self.N);
      for i = 1:self.M
          for j = 1:self.N
              T(i,j,j) = cos(x(j));
              if j == i
                  T(i,j,j) = T(i,j,j) + i*cos(x(i)) + sin(x(i));
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
