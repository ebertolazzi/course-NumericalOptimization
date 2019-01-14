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
        n = int32(2);
      elseif nargin == 1
        n = varargin{1};
      else
        error('Trigonometric: too many arguments in constructor');
      end
      if ~isinteger(n)
        error('Trigonometric: argument must be an integer, found %s',class(n));
      end
      if n <= 1
        error('Trigonometric: argument must be an integer > 1, found %d',n);
      end
      m = n;
      self@FunctionMap(int32(m),int32(n)); 
      self.exact_solutions = zeros(m,0); % unknown solutions, f=0;
      self.guesses         = 1/n*ones(1,n);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      % evaluate the entries (not squared).
      self.check_x(x);
      F    = zeros(m,1);
      for i = 1:m
          F(i) = n + i*(1-cos(x(i))) - sin(x(i));
          for j=1:n
             F(i) = F(i) - cos(x(j)); 
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
      T    = zeros(m,n,n);
      for i = 1:m
          for j = 1:n
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
